-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/mesh_generation/archive/.*$"),
    "Expected this archived script below Model/mesh_generation/archive: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--[[
===============================================================================
SINGLE-NEPHRON SWC -> PROJECTOR-BACKED UGX PIPELINE
===============================================================================

Purpose
-------
Convert one unbranched nephron centerline stored as an SWC file into a coarse
three-domain UG4 mesh:

    Lumen  |  Membrane  |  Inter

The physical interfaces are Apical (Lumen/Membrane), Basolateral
(Membrane/Inter), Inlet/Outlet (the two Lumen caps), and OuterWall (outside of
the padded box). Cap rim edges remain Apical so projector parameters propagate
correctly through repeated refinement.

End-to-end data flow
--------------------
1. Read command-line parameters and derive the output name from the input SWC.
2. Run quadratic_interpolate_swc.py in cubic mode. This softens sharp control
   points and densely resamples the centerline before meshing.
3. Write a temporary SWC whose radius is the OUTER nephron radius:
       outer radius = lumen radius + membrane thickness
4. Call the C++ neuro_collection importer. The importer creates the nested
   Lumen/Membrane tube, a padded box, and a tetrahedral Inter region.
5. Verify that the UGX contains the serialized NeuriteProjector and its
   npSurfParams vertex attachment. These are required for curved refinement.
6. Report geometry counts, diameter, and output path.

Important design decisions
--------------------------
* This script always writes the coarse mesh. Refinement is intentionally left
  to ProMesh so the user begins with a small mesh and obtains progressively
  rounder tube surfaces through the serialized projector.
* The default anisotropy is 6: axial tube cells are longer than their transverse
  size. Pass -anisotropy to override it.
* The Python smoother output, not the raw SWC, is used for mesh generation.
* The default output is runs/<input-file-stem>.ugx. Pass -out to override it.
* Native UG4/TetGen logging is retained. Script-specific debug noise is hidden.

Primary command
---------------
ugshell -ex mesh_generation/archive/swc_to_neurite_projector_nephron_final_pipeline.lua \
  -swc "/path/to/trace.swc"
===============================================================================
]]

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

-- the script directory is used to resolve relative paths for the default SWC and output UGX
-- the input paramaters with call out to the terminal to override the defaults
local script_dir = MODEL_ROOT

-- Configuration --------------------------------------------------------------
-- Every value can be changed from the terminal except numRefs, which is fixed
-- to zero because this production pipeline is deliberately coarse-only.
local cfg = {}
-- Input SWC coordinates and all length parameters below are specified in µm.
-- The temporary SWC passed to the plugin is converted to SI metres.
local MICROMETER_TO_METER = 1.0e-6
cfg.swc = util.GetParam("-swc", script_dir .. "/Nephron traces/3rd_trace_243.swc")
local input_name = cfg.swc:match("([^/]+)$") or cfg.swc
local input_stem = input_name:gsub("%.[^%.]+$", "")
cfg.out = util.GetParam("-out", script_dir .. "/runs/" .. input_stem .. ".ugx")
cfg.smoother = util.GetParam("-smoother",
	script_dir .. "/quadratic_interpolate_swc.py")
cfg.readySwc = util.GetParam("-readySwc", "")
cfg.smoothOrder = math.floor(util.GetParamNumber("-smoothOrder", 3))
cfg.samplesPerCorner = math.floor(util.GetParamNumber("-samplesPerCorner", 30))
cfg.preSmoothAngle = util.GetParamNumber("-preSmoothAngle", 25.0)
cfg.preSmoothStrength = util.GetParamNumber("-preSmoothStrength", 0.85)
cfg.preSmoothPasses = math.floor(util.GetParamNumber("-preSmoothPasses", 35))
cfg.resampleSpacing = util.GetParamNumber("-resampleSpacing", 0.5)
cfg.radius = util.GetParamNumber("-radius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
cfg.anisotropy = util.GetParamNumber("-anisotropy", 6.0)
-- Always save the coarse projector-backed mesh. It can be refined later in ProMesh.
cfg.numRefs = 0
cfg.padding = util.GetParamNumber("-padding", 10.0)
cfg.interQuality = util.GetParamNumber("-interQuality", 2.0)


-- Path, validation, and reporting helpers ------------------------------------
local function dirname(path)
	local dir = path:match("^(.*)/[^/]*$")
	if dir == nil or dir == "" then
		return "."
	end
	return dir
end

local function strip_extension(path)
	return (path:gsub("%.[^/%.]+$", ""))
end

local function shell_quote(s)
	return "'" .. tostring(s):gsub("'", "'\\''") .. "'"
end

local function ensure_dir(path)
	if path ~= "." then
		local ok = os.execute("mkdir -p " .. shell_quote(path))
		assert(ok == true or ok == 0, "Could not create directory: " .. path)
	end
end

local function assert_positive(name, value)
	assert(value > 0.0, name .. " must be > 0")
end

local function radius_tag(radius)
	local tag = string.format("%.8g", radius)
	tag = tag:gsub("%.", "p"):gsub("%-", "m"):gsub("%+", "")
	return tag
end

local function write_radius_override_swc(src, dst, radius, length_scale)
	-- Preserve centerline coordinates and parent connectivity while replacing
	-- every SWC radius with the requested constant outer radius. Root/type
	-- normalization keeps the file compatible with the neurite importer.
	local in_file = assert(io.open(src, "r"), "Cannot open SWC input: " .. src)
	local out_file = assert(io.open(dst, "w"), "Cannot write SWC radius override: " .. dst)

	for line in in_file:lines() do
		local hash_pos = line:find("#", 1, true)
		local payload = line
		local comment = ""
		if hash_pos then
			payload = line:sub(1, hash_pos - 1)
			comment = line:sub(hash_pos)
		end

		local tokens = {}
		for tok in payload:gmatch("%S+") do
			tokens[#tokens + 1] = tok
		end

		if #tokens == 0 then
			out_file:write(line, "\n")
		else
			assert(#tokens == 7, "SWC line must have 7 columns: " .. line)
			local out_line = table.concat({
				tokens[1],
				(tokens[7] == "-1") and "1" or (tokens[2] == "0" and "3" or tokens[2]),
				string.format("%.17g", assert(tonumber(tokens[3])) * length_scale),
				string.format("%.17g", assert(tonumber(tokens[4])) * length_scale),
				string.format("%.17g", assert(tonumber(tokens[5])) * length_scale),
				string.format("%.17g", radius * length_scale),
				tokens[7]
			}, " ")
			if comment ~= "" then
				out_line = out_line .. " " .. comment
			end
			out_file:write(out_line, "\n")
		end
	end

	in_file:close()
	out_file:close()
end

local function file_contains(path, needle)
	local file = assert(io.open(path, "r"), "Cannot open file: " .. path)
	local text = file:read("*a")
	file:close()
	return text:find(needle, 1, true) ~= nil
end

local function count_numbers(text)
	local count = 0
	for _ in text:gmatch("%S+") do count = count + 1 end
	return count
end

local function count_mesh_elements(path)
	-- UGX subset handlers also contain tags named vertices/edges/faces/volumes.
	-- Restrict counting to the connectivity section before <subset_handler>.
	local file = assert(io.open(path, "r"), "Cannot open output UGX: " .. path)
	local text = file:read("*a")
	file:close()
	-- Count only the grid connectivity section. Later subset-handler tags also
	-- use names such as <vertices>, <edges>, <faces>, and <volumes>.
	local geometry = text:match("^(.-)<subset_handler") or text

	local function count_tag(tag, width)
		local total = 0
		for payload in geometry:gmatch("<" .. tag .. "[^>]*>(.-)</" .. tag .. ">") do
			total = total + count_numbers(payload) / width
		end
		return total
	end

	local vertices = count_tag("vertices", 3)
	local edges = count_tag("edges", 2)
	local faces = count_tag("triangles", 3) + count_tag("quadrilaterals", 4)
	local volumes = count_tag("tetrahedrons", 4)
		+ count_tag("pyramids", 5)
		+ count_tag("prisms", 6)
		+ count_tag("hexahedrons", 8)

	return math.floor(vertices + 0.5),
		math.floor(edges + 0.5),
		math.floor(faces + 0.5),
		math.floor(volumes + 0.5)
end

local function assert_neurite_projector_ugx(path)
	assert(file_contains(path, "<projector type=\"neurite\""),
		"UGX is missing a serialized NeuriteProjector: " .. path)
	assert(file_contains(path, "<vertex_attachment name=\"npSurfParams\""),
		"UGX is missing npSurfParams. Refinement cannot project neurite vertices: " .. path)
end

local function refine_with_promesh(src, dst, refs)
	-- Retained as a utility for experiments. The production path below always
	-- uses cfg.numRefs = 0 and therefore does not call this function.
	local mesh = Mesh()
	assert(LoadMesh(mesh, src), "LoadMesh failed: " .. src)
	for i = 1, refs do
		SelectAll(mesh)
		print("  projected refine " .. i .. " / " .. refs)
		Refine(mesh)
	end
	assert(SaveMesh(mesh, dst), "SaveMesh failed: " .. dst)
end

-- Validate all user-controlled geometry and smoothing parameters -------------
assert(cfg.swc ~= "", "Pass an SWC file with -swc <path>")
assert(cfg.out ~= "", "Pass an output path with -out <path>")
assert_positive("anisotropy", cfg.anisotropy)
assert_positive("padding", cfg.padding)
assert_positive("interQuality", cfg.interQuality)
assert(cfg.smoothOrder >= 1 and cfg.smoothOrder <= 3, "smoothOrder must be 1, 2, or 3")
assert(cfg.samplesPerCorner >= 1, "samplesPerCorner must be >= 1")
assert(cfg.preSmoothAngle >= 0.0, "preSmoothAngle must be >= 0")
assert(cfg.preSmoothStrength >= 0.0 and cfg.preSmoothStrength <= 1.0,
	"preSmoothStrength must be between 0 and 1")
assert(cfg.preSmoothPasses >= 0, "preSmoothPasses must be >= 0")
assert(cfg.resampleSpacing >= 0.0, "resampleSpacing must be >= 0")
assert_positive("membraneThickness", cfg.membraneThickness)
assert(cfg.numRefs >= 0, "numRefs must be >= 0")
if cfg.radius == 0.0 then
	error("radius must be > 0, or < 0 to keep the SWC radii")
end

ensure_dir(dirname(cfg.out))

-- Stage 1: smooth and resample the raw centerline -----------------------------
local ready_swc = cfg.readySwc
if ready_swc == "" then
	ready_swc = strip_extension(cfg.swc) .. "_mesh_ready.swc"
end
ensure_dir(dirname(ready_swc))

local smooth_command = table.concat({
	"python3",
	shell_quote(cfg.smoother),
	shell_quote(cfg.swc),
	"-o", shell_quote(ready_swc),
	"--order", tostring(cfg.smoothOrder),
	"--samples-per-corner", tostring(cfg.samplesPerCorner),
	"--pre-smooth-angle", tostring(cfg.preSmoothAngle),
	"--pre-smooth-strength", tostring(cfg.preSmoothStrength),
	"--pre-smooth-passes", tostring(cfg.preSmoothPasses),
	"--resample-spacing", tostring(cfg.resampleSpacing),
	"> /dev/null"
}, " ")

local smooth_ok = os.execute(smooth_command)
assert(smooth_ok == true or smooth_ok == 0,
	"quadratic_interpolate_swc.py failed")

-- Stage 2: encode the constant outer radius in a temporary SWC ----------------
-- The C++ nested-tube generator treats the SWC radius as the outer tube radius.
-- lumen_scale below then places the inner Lumen/Apical surface at radius 1.
local import_swc = ready_swc
local outer_radius = cfg.radius
if cfg.radius > 0.0 then
	assert_positive("radius", cfg.radius)
	outer_radius = cfg.radius + cfg.membraneThickness
	import_swc = strip_extension(cfg.out) .. "_outer_radius_" .. radius_tag(outer_radius) .. ".swc"
	-- Build the mesh from the interpolated/smoothed trace.  Using cfg.swc here
	-- silently discarded the smoothing stage and could feed sharp raw turns to
	-- TetGen.
	write_radius_override_swc(
		ready_swc, import_swc, outer_radius, MICROMETER_TO_METER)
else
	error("The membrane test currently requires -radius > 0 so it can apply a constant physical membrane thickness.")
end
local lumen_scale = cfg.radius / outer_radius

local out_ugx = strip_extension(cfg.out) .. ".ugx"
local coarse_ugx = out_ugx

-- Stage 3: create Lumen, Membrane, Inter, interfaces, box, and projector ------
-- Keep importer refinement disabled here. Its hierarchy files are useful for
-- some solvers, but this workflow saves one clean coarse UGX with defPH.
import_nephron_with_membrane_and_inter_from_swc(import_swc, coarse_ugx,
	lumen_scale, cfg.anisotropy, 0,
	cfg.padding * MICROMETER_TO_METER, cfg.interQuality)
assert_neurite_projector_ugx(coarse_ugx)

-- Stage 4: concise result summary --------------------------------------------
local numVertices, numEdges, numFaces, numVolumes = count_mesh_elements(out_ugx)
print("Vertices: " .. numVertices)
print("Edges: " .. numEdges)
print("Faces: " .. numFaces)
print("Volumes: " .. numVolumes)
print("Diameter (m): " .. (2.0 * cfg.radius * MICROMETER_TO_METER))
print("Output: " .. out_ugx)
