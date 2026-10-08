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
MULTIPLE-NEPHRON SWC -> ONE SHARED PROJECTOR-BACKED UGX PIPELINE
===============================================================================

Purpose
-------
Build two or more disconnected nephron traces inside one conforming padded box.
Each nephron keeps independent biological subsets, while all traces share one
extracellular/interstitial domain:

Per input i:
    Membrane_i, Lumen_i, Basolateral_i, Apical_i, Inlet_i, Outlet_i

Shared:
    Inter, OuterWall

The input order controls suffix numbering: -swc1 creates *_1, -swc2 creates
*_2, and so on.

End-to-end data flow
--------------------
1. Collect the contiguous command-line sequence -swc1, -swc2, ... .
2. Smooth every input independently with quadratic_interpolate_swc.py.
3. Read the smoothed traces, replace their radii with the constant outer radius,
   and remap node/parent IDs into one collision-free, multi-root SWC.
4. Call the dedicated C++ multi-nephron importer exactly once. Calling the
   single-nephron importer repeatedly would create overlapping boxes and
   nonconforming Inter domains.
5. The C++ importer generates all disconnected nested tubes, one shared box,
   one shared TetGen Inter mesh, indexed subsets, terminal caps, and projector
   mappings for every Apical_i/Basolateral_i surface.
6. Verify projector serialization and print mesh counts.

Topology and refinement guarantees
----------------------------------
* Each input SWC must be one unbranched root-to-terminal path.
* Each root end becomes Inlet_i; the opposite axial end becomes Outlet_i.
* Each cap face is Inlet_i/Outlet_i. Its shared rim edges remain Apical_i so
  NeuriteProjector can propagate npSurfParams through repeated refinement.
* Only Apical_i and Basolateral_i receive NeuriteProjectors. Volumes, caps,
  Inter, and OuterWall are not projected onto a tube.
* The output is coarse-only and can be refined later in ProMesh.

Primary command
---------------
ugshell -ex mesh_generation/archive/swc_to_neurite_projector_multiple_nephrons.lua \
  -swc1 "/path/first.swc" -swc2 "/path/second.swc" \
  -out "./runs/multiple_nephrons.ugx"
===============================================================================
]]

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

-- Configuration --------------------------------------------------------------
local script_dir = MODEL_ROOT
local cfg = {}
-- Input SWC coordinates and all length parameters are in µm. The merged SWC
-- passed to the plugin is converted to SI metres.
local MICROMETER_TO_METER = 1.0e-6
cfg.out = util.GetParam("-out", script_dir .. "/runs/multiple_nephrons.ugx")
cfg.smoother = util.GetParam("-smoother",
	script_dir .. "/quadratic_interpolate_swc.py")
cfg.radius = util.GetParamNumber("-radius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
cfg.anisotropy = util.GetParamNumber("-anisotropy", 6.0)
cfg.padding = util.GetParamNumber("-padding", 10.0)
cfg.interQuality = util.GetParamNumber("-interQuality", 2.0)
cfg.smoothOrder = math.floor(util.GetParamNumber("-smoothOrder", 3))
cfg.samplesPerCorner = math.floor(util.GetParamNumber("-samplesPerCorner", 30))
cfg.preSmoothAngle = util.GetParamNumber("-preSmoothAngle", 25.0)
cfg.preSmoothStrength = util.GetParamNumber("-preSmoothStrength", 0.85)
cfg.preSmoothPasses = math.floor(util.GetParamNumber("-preSmoothPasses", 35))
cfg.resampleSpacing = util.GetParamNumber("-resampleSpacing", 0.5)

-- Input discovery ------------------------------------------------------------
-- Stop at the first missing index. Do not skip numbers (for example, do not
-- provide -swc1 and -swc3 without -swc2).
local swcs = {}
for i = 1, 1000 do
	local path = util.GetParam("-swc" .. i, "")
	if path == "" then break end
	swcs[#swcs + 1] = path
end
assert(#swcs >= 2, "Pass at least -swc1 <file> and -swc2 <file>.")

local function dirname(path)
	return path:match("^(.*)/[^/]*$") or "."
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

local function count_numbers(text)
	local count = 0
	for _ in text:gmatch("%S+") do count = count + 1 end
	return count
end

local function count_mesh_elements(path)
	-- Count grid connectivity only; subset index lists reuse the same XML names.
	local file = assert(io.open(path, "r"), "Cannot open output UGX: " .. path)
	local text = file:read("*a")
	file:close()
	local geometry = text:match("^(.-)<subset_handler") or text
	local function count_tag(tag, width)
		local total = 0
		for payload in geometry:gmatch("<" .. tag .. "[^>]*>(.-)</" .. tag .. ">") do
			total = total + count_numbers(payload) / width
		end
		return total
	end
	return math.floor(count_tag("vertices", 3) + 0.5),
		math.floor(count_tag("edges", 2) + 0.5),
		math.floor(count_tag("triangles", 3) + count_tag("quadrilaterals", 4) + 0.5),
		math.floor(count_tag("tetrahedrons", 4) + count_tag("pyramids", 5)
			+ count_tag("prisms", 6) + count_tag("hexahedrons", 8) + 0.5)
end

local function file_contains(path, needle)
	local file = assert(io.open(path, "r"), "Cannot open file: " .. path)
	local text = file:read("*a")
	file:close()
	return text:find(needle, 1, true) ~= nil
end

local function smooth_swc(src, dst)
	-- Smooth traces independently so interpolation never bridges two nephrons.
	local command = table.concat({
		"python3", shell_quote(cfg.smoother), shell_quote(src),
		"-o", shell_quote(dst),
		"--order", tostring(cfg.smoothOrder),
		"--samples-per-corner", tostring(cfg.samplesPerCorner),
		"--pre-smooth-angle", tostring(cfg.preSmoothAngle),
		"--pre-smooth-strength", tostring(cfg.preSmoothStrength),
		"--pre-smooth-passes", tostring(cfg.preSmoothPasses),
		"--resample-spacing", tostring(cfg.resampleSpacing),
		"> /dev/null"
	}, " ")
	local ok = os.execute(command)
	assert(ok == true or ok == 0, "SWC smoothing failed: " .. src)
end

local function read_swc(path)
	-- Parse the seven standard SWC fields. Parent consistency is checked later
	-- while IDs are remapped into the merged multi-root file.
	local points = {}
	local file = assert(io.open(path, "r"), "Cannot open SWC: " .. path)
	for line in file:lines() do
		local payload = line:gsub("#.*$", "")
		local tokens = {}
		for token in payload:gmatch("%S+") do tokens[#tokens + 1] = token end
		if #tokens > 0 then
			assert(#tokens == 7, "SWC line must have 7 columns: " .. line)
			points[#points + 1] = tokens
		end
	end
	file:close()
	return points
end

assert(cfg.radius > 0, "radius must be > 0")
assert(cfg.membraneThickness > 0, "membraneThickness must be > 0")
assert(cfg.anisotropy > 0, "anisotropy must be > 0")
assert(cfg.padding > 0, "padding must be > 0")
assert(cfg.interQuality > 0, "interQuality must be > 0")
ensure_dir(dirname(cfg.out))

-- Stage 1: smooth, radius-normalize, and merge all traces ---------------------
local out_base = strip_extension(cfg.out)
local merged_swc = out_base .. "_merged_outer_radius.swc"
local merged = assert(io.open(merged_swc, "w"), "Cannot create merged SWC.")
local next_id = 1
local outer_radius = cfg.radius + cfg.membraneThickness

for trace_index, raw_swc in ipairs(swcs) do
	local ready_swc = strip_extension(raw_swc) .. "_mesh_ready.swc"
	ensure_dir(dirname(ready_swc))
	smooth_swc(raw_swc, ready_swc)
	local points = read_swc(ready_swc)
	-- Every source file may reuse IDs such as 1, 2, 3. Allocate globally unique
	-- IDs and translate each non-root parent through this file-local map.
	local id_map = {}
	for _, point in ipairs(points) do
		id_map[tonumber(point[1])] = next_id
		next_id = next_id + 1
	end
	for _, point in ipairs(points) do
		local old_parent = tonumber(point[7])
		local new_parent = old_parent == -1 and -1 or assert(id_map[old_parent],
			"Missing parent in " .. ready_swc)
		local point_type = new_parent == -1 and 1
			or (tonumber(point[2]) == 0 and 3 or tonumber(point[2]))
		merged:write(table.concat({
			id_map[tonumber(point[1])], point_type,
			string.format("%.17g", assert(tonumber(point[3])) * MICROMETER_TO_METER),
			string.format("%.17g", assert(tonumber(point[4])) * MICROMETER_TO_METER),
			string.format("%.17g", assert(tonumber(point[5])) * MICROMETER_TO_METER),
			string.format("%.17g", outer_radius * MICROMETER_TO_METER), new_parent
		}, " "), "\n")
	end
end
merged:close()

-- Stage 2: create all nested tubes plus one shared Inter/OuterWall domain -----
local out_ugx = out_base .. ".ugx"
local lumen_scale = cfg.radius / outer_radius
import_multiple_nephrons_with_membrane_and_inter_from_swc(
	merged_swc, out_ugx, #swcs, lumen_scale, cfg.anisotropy,
	0, cfg.padding * MICROMETER_TO_METER, cfg.interQuality)

-- Stage 3: fail early if refinement metadata was not serialized ---------------
assert(file_contains(out_ugx, "<projector type=\"neurite\""),
	"UGX is missing its NeuriteProjector.")
assert(file_contains(out_ugx, "<vertex_attachment name=\"npSurfParams\""),
	"UGX is missing npSurfParams.")

-- Stage 4: concise result summary --------------------------------------------
local vertices, edges, faces, volumes = count_mesh_elements(out_ugx)
print("Nephrons: " .. #swcs)
print("Vertices: " .. vertices)
print("Edges: " .. edges)
print("Faces: " .. faces)
print("Volumes: " .. volumes)
print("Diameter (m): " .. (2.0 * cfg.radius * MICROMETER_TO_METER))
print("Output: " .. out_ugx)
