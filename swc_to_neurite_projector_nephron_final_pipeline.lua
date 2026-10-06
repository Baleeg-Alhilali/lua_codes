-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Build a projector-backed neurite surface directly from the specified SWC tree.
--
-- This is a copy of swc_to_neurite_projector.lua tailored to the
-- Nephron traces dataset and a fixed radius of 1.0.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

-- the script directory is used to resolve relative paths for the default SWC and output UGX
-- the input paramaters with call out to the terminal to override the defaults
local script_dir = MODEL_ROOT

local cfg = {}
-- SWC coordinates and CLI lengths are in um; UGX/projector coordinates are mm.
local MESH_LENGTH_SCALE = 1e-3
cfg.swc = util.GetParam("-swc", script_dir .. "/Nephron traces/3rd_trace_243.swc")
local input_name = cfg.swc:match("([^/]+)$") or cfg.swc
local input_stem = input_name:gsub("%.[^%.]+$", "")
cfg.out = util.GetParam("-out", script_dir .. "/runs/" .. input_stem .. ".ugx")
cfg.smoother = util.GetParam("-smoother",
	script_dir .. "/quadratic_interpolate_swc.py")
cfg.readySwc = util.GetParam("-readySwc", "")
cfg.extendTerminals = util.GetParamNumber("-extendTerminals", 1) ~= 0
cfg.terminalExtender = util.GetParam("-terminalExtender",
	script_dir .. "/extend_swc_terminals_to_outerwall.py")
cfg.terminalPaddingFactor = util.GetParamNumber("-terminalPaddingFactor", 25.0)
cfg.smoothOrder = math.floor(util.GetParamNumber("-smoothOrder", 3))
cfg.samplesPerCorner = math.floor(util.GetParamNumber("-samplesPerCorner", 30))
-- Cubic interpolation can reintroduce 20--25 degree turns even after the raw
-- controls were softened at 25 degrees.  A 10-degree control threshold keeps
-- the coarse axial sweep valid when refinement uses the exact projector.
cfg.preSmoothAngle = util.GetParamNumber("-preSmoothAngle", 10.0)
cfg.preSmoothStrength = util.GetParamNumber("-preSmoothStrength", 0.85)
cfg.preSmoothPasses = math.floor(util.GetParamNumber("-preSmoothPasses", 50))
cfg.resampleSpacing = util.GetParamNumber("-resampleSpacing", 0.5)
cfg.radius = util.GetParamNumber("-radius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
cfg.anisotropy = util.GetParamNumber("-anisotropy", 1.0)
cfg.ogrid = math.floor(util.GetParamNumber("-ogrid", 10))
-- Keep one lumen-midpoint support ring between the center and Apical ring.
-- Removing it makes projected cross-sectional refinement connect large child
-- triangles directly to alternating coarse angular vertices.  The reduced
-- center-only topology remains available with -coarseLumenCenterOnly 1.
cfg.coarseLumenCenterOnly = util.GetParamNumber("-coarseLumenCenterOnly", 0) ~= 0
-- Refine the complete domain by default.  -numRefs N writes standalone flat
-- inspection meshes _refined_1.ugx through _refined_N.ugx.  The flow solver
-- loads the coarse projector-backed UGX and builds its multigrid hierarchy in
-- memory. Cross-sectional-only refinement remains experimental (-crossOnly 1).
cfg.numRefs = math.floor(util.GetParamNumber("-numRefs", 0))
cfg.crossOnly = util.GetParamNumber("-crossOnly", 0) ~= 0
cfg.padding = util.GetParamNumber("-padding", 10.0)
cfg.interQuality = util.GetParamNumber("-interQuality", 2.0)
-- Prefer -BoxRefs, while retaining the older -boxRefs spelling as an alias.
cfg.boxRefs = math.floor(util.GetParamNumber("-BoxRefs",
	util.GetParamNumber("-boxRefs", 4)))


-- Utility functions
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

local function write_radius_override_swc(src, dst, radius)
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
				string.format("%.17g", assert(tonumber(tokens[3])) * MESH_LENGTH_SCALE),
				string.format("%.17g", assert(tonumber(tokens[4])) * MESH_LENGTH_SCALE),
				string.format("%.17g", assert(tonumber(tokens[5])) * MESH_LENGTH_SCALE),
				string.format("%.17g", radius * MESH_LENGTH_SCALE),
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
	assert(file_contains(path, "<subset name=\"TripleJunction\""),
		"UGX is missing the terminal TripleJunction subset: " .. path)
end

local function refine_with_promesh(src, dst, refs)
	local mesh = Mesh()
	assert(LoadMesh(mesh, src), "LoadMesh failed: " .. src)
	for i = 1, refs do
		SelectAll(mesh)
		print("  projected refine " .. i .. " / " .. refs)
		Refine(mesh)
	end
	assert(SaveMesh(mesh, dst), "SaveMesh failed: " .. dst)
end

assert(cfg.swc ~= "", "Pass an SWC file with -swc <path>")
assert(cfg.out ~= "", "Pass an output path with -out <path>")
assert_positive("anisotropy", cfg.anisotropy)
assert(math.abs(cfg.anisotropy - 1.0) < 1e-12,
	"anisotropy must be 1 for the verified projector-safe nephron pipeline")
assert(cfg.ogrid >= 4, "ogrid must be an integer >= 4")
assert_positive("padding", cfg.padding)
assert(cfg.terminalPaddingFactor >= 0.0, "terminalPaddingFactor must be >= 0")
assert_positive("interQuality", cfg.interQuality)
assert(cfg.boxRefs >= 0 and cfg.boxRefs <= 10,
	"BoxRefs must be an integer from 0 to 10")
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

local import_swc = ready_swc
local outer_radius = cfg.radius
if cfg.radius > 0.0 then
	assert_positive("radius", cfg.radius)
	outer_radius = cfg.radius + cfg.membraneThickness
	import_swc = strip_extension(cfg.out) .. "_outer_radius_" .. radius_tag(outer_radius) .. ".swc"
	-- Build the mesh from the interpolated/smoothed trace.  Using cfg.swc here
	-- silently discarded the smoothing stage and could feed sharp raw turns to
	-- TetGen.
	write_radius_override_swc(ready_swc, import_swc, outer_radius)
else
	error("The membrane test currently requires -radius > 0 so it can apply a constant physical membrane thickness.")
end
local lumen_scale = cfg.radius / outer_radius

if cfg.extendTerminals then
	local extended_swc = strip_extension(cfg.out) .. "_outerwall_extended.swc"
	local extension_command = table.concat({
		"python3", shell_quote(cfg.terminalExtender), shell_quote(import_swc),
		"-o", shell_quote(extended_swc),
		"--padding", tostring(cfg.padding * MESH_LENGTH_SCALE),
		"--spacing", tostring(math.max(cfg.resampleSpacing, 0.1) * MESH_LENGTH_SCALE),
		"--box-refs", tostring(cfg.boxRefs),
		"--min-padding-factor", tostring(cfg.terminalPaddingFactor)
	}, " ")
	local extension_ok = os.execute(extension_command)
	assert(extension_ok == true or extension_ok == 0,
		"Terminal-to-OuterWall SWC extension failed")
	import_swc = extended_swc
end

local out_ugx = strip_extension(cfg.out) .. ".ugx"
local coarse_ugx = out_ugx

-- Build the standard single-nephron Lumen/Membrane/Inter geometry. Terminal
-- cap-rim edges and vertices are placed in TripleJunction so the three
-- compartment fields have nodal support without extending Inter onto Apical.
-- The
-- importer performs each requested refinement with its serialized neurite
-- projector and writes a separate single-level UGX for visual inspection.
import_nephron_with_membrane_and_inter_from_swc(import_swc, coarse_ugx,
	lumen_scale, cfg.anisotropy, cfg.numRefs, cfg.padding * MESH_LENGTH_SCALE,
	cfg.interQuality, cfg.boxRefs, cfg.ogrid, cfg.coarseLumenCenterOnly,
	cfg.crossOnly)
assert_neurite_projector_ugx(coarse_ugx)

local numVertices, numEdges, numFaces, numVolumes = count_mesh_elements(out_ugx)
print("Vertices: " .. numVertices)
print("Edges: " .. numEdges)
print("Faces: " .. numFaces)
print("Volumes: " .. numVolumes)
print("Lumen diameter (um): " .. (2.0 * cfg.radius))
print("Mesh coordinate unit: mm")
print("Terminal-to-OuterWall extension: " .. tostring(cfg.extendTerminals))
print("Output: " .. out_ugx)
