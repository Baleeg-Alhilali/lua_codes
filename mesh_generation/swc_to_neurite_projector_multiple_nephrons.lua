-- Resolve inputs and outputs from Model while allowing this script to be run
-- directly from Model/mesh_generation.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
	__model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local GENERATION_DIR = assert(__model_script:match("^(.*)/[^/]+$"),
	"Cannot determine mesh_generation directory: " .. __model_script)
local MODEL_ROOT = assert(GENERATION_DIR:match("^(.*)/mesh_generation$"),
	"Expected this script inside Model/mesh_generation: " .. __model_script)
ChangeDirectory(GENERATION_DIR)

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
ugshell -ex swc_to_neurite_projector_multiple_nephrons.lua \
  -swc1 "/path/first.swc" -swc2 "/path/second.swc" \
  -out "./runs/multiple_nephrons.ugx"

True-scale workflow (2026-11-24, 1 um voxels; see
docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md and
tools/rod_repair_joint_workflow.sh):
1. repair the traces jointly, length-preserving, with inter-trace clearance:
   tools/repair_trace_rod.py 3rd.swc 4th.swc --free-ends -o a.swc b.swc
2. this script with -swc1 a.swc -swc2 b.swc -skipSmoothing 1 -numRefs 0
   writes <out>_outerwall_extended.swc (merged, multi-root, mm);
3. repair that file jointly with frozen OuterWall tails:
   tools/repair_trace_rod.py <out>_outerwall_extended.swc --voxel-size 1000
       --frozen-head 40 --frozen-tail 40 -o merged_repaired.swc
4. final mesh: -mergedSwc merged_repaired.swc -interMesh hex
   -prismInterThickness 3 -ogrid 12 -BoxRefs 3 -numRefs 2
   (no smoothing, no further extension; the UG4_ box metadata is reused).
===============================================================================
]]

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

-- Configuration --------------------------------------------------------------
local script_dir = GENERATION_DIR
local cfg = {}
-- Input SWC coordinates and all CLI lengths are in um; UGX/projector
-- coordinates are in mm. Smoothing is done before this conversion.
local MESH_LENGTH_SCALE = 1e-3
cfg.out = util.GetParam("-out", MODEL_ROOT .. "/mesh_generation/meshes/" .. os.date("%Y-%m-%d") .. "_multiple_nephrons/multiple_nephrons.ugx")
cfg.smoother = util.GetParam("-smoother",
	script_dir .. "/quadratic_interpolate_swc.py")
-- Copy every -swcN trace verbatim instead of smoothing it (use for traces
-- already repaired by tools/repair_trace_rod.py, in micrometres).
cfg.skipSmoothing = util.GetParamNumber("-skipSmoothing", 0) ~= 0
-- A pre-merged multi-root SWC in micrometres (for example the joint rod
-- repair of this script's own _outerwall_extended.swc).  Replaces -swcN:
-- no smoothing, no merging and no terminal extension; its # UG4_ box
-- metadata is carried into the import file.
cfg.mergedSwc = util.GetParam("-mergedSwc", "")
cfg.extendTerminals = util.GetParamNumber("-extendTerminals", 1) ~= 0
cfg.terminalExtender = util.GetParam("-terminalExtender",
	script_dir .. "/extend_swc_terminals_to_outerwall.py")
-- Same defaults as the single-nephron pipeline (2026-11-24): 25 x outer radius
-- made a 500 um minimum box padding at true scale.
cfg.terminalPaddingFactor = util.GetParamNumber("-terminalPaddingFactor", 0.0)
-- True-scale defaults (1 um voxels): 10 um lumen radius, 10 um membrane.
cfg.radius = util.GetParamNumber("-radius", 10.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 10.0)
cfg.centerlineScale = util.GetParamNumber("-centerlineScale", 1.0)
local requested_outer_radius = cfg.radius + cfg.membraneThickness
cfg.skipGeometryPreflight = util.GetParamNumber("-skipGeometryPreflight", 0) ~= 0
cfg.anisotropy = util.GetParamNumber("-anisotropy", 4.0)
cfg.adaptiveAnisotropy = util.GetParamNumber("-adaptiveAnisotropy", 1) ~= 0
cfg.bendAnisotropy = util.GetParamNumber("-bendAnisotropy", 1.0)
cfg.bendChordErrorFactor = util.GetParamNumber("-bendChordErrorFactor", 0.01)
cfg.gradedInter = util.GetParamNumber("-gradedInter", 1) ~= 0
cfg.interNearEdgeFactor = util.GetParamNumber("-interNearEdgeFactor", 2.0)
cfg.interNearDistanceFactor = util.GetParamNumber("-interNearDistanceFactor", 1.5)
cfg.interFarDistanceFactor = util.GetParamNumber("-interFarDistanceFactor", 6.0)
cfg.interFarEdgeFactor = util.GetParamNumber("-interFarEdgeFactor", 6.0)
cfg.ogrid = math.floor(util.GetParamNumber("-ogrid", 10))
-- Near-Inter topology, identical to the single-nephron pipeline: "hex" builds
-- an all-hex structured shell of -prismInterThickness um with outward pyramid
-- caps (which refine linearly; UG_NEPHRON_PROJECTED_CAPS=1 restores projected
-- caps); "prism" keeps the older prismatic shell.
cfg.interMesh = util.GetParam("-interMesh", "prism")
cfg.prismInterThickness = util.GetParamNumber("-prismInterThickness", 20.0)
cfg.interFadeDistance = util.GetParamNumber("-interFadeDistance", 0.0)
cfg.coordinateRelTol = util.GetParamNumber("-coordinateRelTol", 1e-10)
cfg.containmentRelTol = util.GetParamNumber("-containmentRelTol", 1e-6)
cfg.parameterTol = util.GetParamNumber("-parameterTol", 1e-7)
cfg.mergeRadiusFactor = util.GetParamNumber("-mergeRadiusFactor", 1e-4)
cfg.terminalSnapRadiusFactor = util.GetParamNumber("-terminalSnapRadiusFactor", 2.0)
cfg.straightnessRadiusFactor = util.GetParamNumber("-straightnessRadiusFactor", 1e-7)
cfg.angularSpacingToleranceFactor = util.GetParamNumber(
	"-angularSpacingToleranceFactor", 1e-5)
cfg.directionTolerance = util.GetParamNumber("-directionTolerance", 1e-8)
cfg.interMinLayers = math.floor(util.GetParamNumber("-interMinLayers", 4))
cfg.interMaxLayers = math.floor(util.GetParamNumber("-interMaxLayers", 12))
cfg.interLayersPerSurfaceEdge = util.GetParamNumber("-interLayersPerSurfaceEdge", 3.0)
cfg.interRadialGrowth = util.GetParamNumber("-interRadialGrowth",
	cfg.interMesh == "hex" and 0.0 or 2.0)
-- Structural validation of the coarse and every refined file, as in the
-- single-nephron pipeline.
cfg.validate = util.GetParamNumber("-validate", 1) ~= 0
cfg.validator = util.GetParam("-validator", script_dir .. "/validate_ugx_mesh.py")
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
cfg.padding = util.GetParamNumber("-padding",
	math.max(10.0 * cfg.centerlineScale, 2.0 * requested_outer_radius))
cfg.paddingX = util.GetParamNumber("-paddingX", cfg.padding)
cfg.paddingY = util.GetParamNumber("-paddingY", cfg.padding)
cfg.paddingZ = util.GetParamNumber("-paddingZ", cfg.padding)
cfg.interQuality = util.GetParamNumber("-interQuality", 2.0)
-- Prefer -BoxRefs, while retaining the older -boxRefs spelling as an alias.
cfg.boxRefs = math.floor(util.GetParamNumber("-BoxRefs",
	util.GetParamNumber("-boxRefs", 4)))
cfg.smoothOrder = math.floor(util.GetParamNumber("-smoothOrder", 3))
cfg.samplesPerCorner = math.floor(util.GetParamNumber("-samplesPerCorner", 30))
-- Match the verified single-nephron coarse-sweep repair.  This prevents cubic
-- interpolation from leaving sharp local turns that fold exact projected hexes.
cfg.preSmoothAngle = util.GetParamNumber("-preSmoothAngle", 10.0)
cfg.preSmoothStrength = util.GetParamNumber("-preSmoothStrength", 0.85)
cfg.preSmoothPasses = math.floor(util.GetParamNumber("-preSmoothPasses", 50))
cfg.resampleSpacing = util.GetParamNumber("-resampleSpacing", 5.0)
cfg.terminalSpacing = util.GetParamNumber("-terminalSpacing", 20.0)
-- Minimum distance (um) between terminal ports of different nephrons on the
-- same OuterWall face; default 3 x outer radius (60 um for 10 + 10 um).
cfg.minPortSeparation = util.GetParamNumber("-minPortSeparation",
	3.0 * requested_outer_radius)

-- Input discovery ------------------------------------------------------------
-- Stop at the first missing index. Do not skip numbers (for example, do not
-- provide -swc1 and -swc3 without -swc2).
local swcs = {}
for i = 1, 1000 do
	local path = util.GetParam("-swc" .. i, "")
	if path == "" then break end
	swcs[#swcs + 1] = path
end
if cfg.mergedSwc ~= "" then
	assert(#swcs == 0, "Pass either -mergedSwc or -swc1 -swc2 ..., not both.")
else
	assert(#swcs >= 2, "Pass at least -swc1 <file> and -swc2 <file>, or -mergedSwc <file>.")
end

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
	if cfg.skipSmoothing then
		local i = assert(io.open(src, "r"), "Cannot open SWC: " .. src)
		local text = i:read("*a")
		i:close()
		local o = assert(io.open(dst, "w"), "Cannot write SWC: " .. dst)
		o:write(text)
		o:close()
		print("Smoothing skipped: using " .. src .. " unchanged")
		return
	end
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

local function audit_trace_curvature(points, label, outer_radius)
	assert(#points >= 3, "Geometry preflight needs at least three points: " .. label)
	local function xyz(point)
		return assert(tonumber(point[3])) * cfg.centerlineScale,
			assert(tonumber(point[4])) * cfg.centerlineScale,
			assert(tonumber(point[5])) * cfg.centerlineScale
	end
	local function distance(a, b)
		local ax, ay, az = xyz(a)
		local bx, by, bz = xyz(b)
		local dx, dy, dz = ax-bx, ay-by, az-bz
		return math.sqrt(dx*dx + dy*dy + dz*dz)
	end
	local minimum = math.huge
	for i = 2, #points-1 do
		local ax, ay, az = xyz(points[i-1])
		local bx, by, bz = xyz(points[i])
		local cx, cy, cz = xyz(points[i+1])
		local abx, aby, abz = bx-ax, by-ay, bz-az
		local acx, acy, acz = cx-ax, cy-ay, cz-az
		local nx = aby*acz-abz*acy
		local ny = abz*acx-abx*acz
		local nz = abx*acy-aby*acx
		local twice_area = math.sqrt(nx*nx + ny*ny + nz*nz)
		if twice_area > 1e-12 * cfg.centerlineScale * cfg.centerlineScale then
			minimum = math.min(minimum,
				distance(points[i-1], points[i])
				* distance(points[i], points[i+1])
				* distance(points[i-1], points[i+1]) / (2.0 * twice_area))
		end
	end
	print(string.format("Geometry preflight %s (um): outerRadius=%.9g, minCurvatureRadius=%.9g",
		label, outer_radius, minimum))
	if not cfg.skipGeometryPreflight then
		assert(outer_radius < 0.9 * minimum,
			string.format("%s: requested outer radius %.6g um is incompatible with centerline curvature %.6g um. Reduce radius/thickness or increase -centerlineScale.",
				label, outer_radius, minimum))
	end
end

local function validate_ugx(path)
	if not cfg.validate then return end
	print("Validating generated mesh: " .. path)
	local ok = os.execute("python3 " .. shell_quote(cfg.validator) .. " " .. shell_quote(path))
	assert(ok == true or ok == 0, "Mesh validation failed: " .. path)
end

-- Non-local clearance within each trace and the distance between traces
-- (um), on the points as given.  Traces closer than 2 x outer / 0.9 overlap
-- or leave no room for the shared Inter mesh.
local function audit_clearance(traces, outer_radius)
	local function d(a, b)
		local dx, dy, dz = a[1]-b[1], a[2]-b[2], a[3]-b[3]
		return math.sqrt(dx*dx + dy*dy + dz*dz)
	end
	local exclusion = math.max(4.0 * outer_radius, 20.0)
	local self_min, cross_min = math.huge, math.huge
	for t, pts in ipairs(traces) do
		local s = {0.0}
		for i = 2, #pts do s[i] = s[i-1] + d(pts[i-1], pts[i]) end
		for i = 1, #pts do
			for j = i + 1, #pts do
				if s[j] - s[i] > exclusion then self_min = math.min(self_min, d(pts[i], pts[j])) end
			end
			for u = t + 1, #traces do
				for _, q in ipairs(traces[u]) do cross_min = math.min(cross_min, d(pts[i], q)) end
			end
		end
	end
	print(string.format("Geometry preflight (um): minNonlocalCenterlineDistance=%.9g, minDistanceBetweenTraces=%.9g",
		self_min, cross_min))
	if not cfg.skipGeometryPreflight then
		assert(2.0 * outer_radius < 0.9 * math.min(self_min, cross_min),
			"Traces (or distant parts of one trace) are closer than 2 x outer radius / 0.9; "
			.. "repair them jointly with tools/repair_trace_rod.py.")
	end
end

local function trace_points(points, scale)
	local pts = {}
	for _, p in ipairs(points) do
		pts[#pts + 1] = {assert(tonumber(p[3])) * scale, assert(tonumber(p[4])) * scale,
			assert(tonumber(p[5])) * scale}
	end
	return pts
end

assert(cfg.interMesh == "hex" or cfg.interMesh == "prism", "-interMesh must be 'hex' or 'prism'")
assert(cfg.prismInterThickness > 0, "prismInterThickness must be > 0")
assert(cfg.interFadeDistance >= 0, "interFadeDistance must be >= 0")
assert(cfg.radius > 0, "radius must be > 0")
assert(cfg.membraneThickness > 0, "membraneThickness must be > 0")
assert(cfg.centerlineScale > 0, "centerlineScale must be > 0")
assert(cfg.terminalSpacing > 0, "terminalSpacing must be > 0")
assert(cfg.anisotropy > 0, "anisotropy must be > 0")
assert(cfg.anisotropy >= 1.0, "anisotropy must be >= 1")
assert(cfg.bendAnisotropy >= 1.0 and cfg.bendAnisotropy <= cfg.anisotropy,
	"bendAnisotropy must be between 1 and anisotropy")
assert(cfg.bendChordErrorFactor > 0.0, "bendChordErrorFactor must be > 0")
assert(cfg.interNearEdgeFactor > 0.0, "interNearEdgeFactor must be > 0")
assert(cfg.interFarEdgeFactor >= cfg.interNearEdgeFactor,
	"interFarEdgeFactor must be >= interNearEdgeFactor")
assert(cfg.interNearDistanceFactor >= 0.0
	and cfg.interFarDistanceFactor > cfg.interNearDistanceFactor,
	"Inter distance factors require 0 <= near < far")
assert(cfg.ogrid >= 4, "ogrid must be an integer >= 4")
assert(cfg.numRefs >= 0, "numRefs must be >= 0")
assert(cfg.padding > 0, "padding must be > 0")
assert(cfg.paddingX > 0, "paddingX must be > 0")
assert(cfg.paddingY > 0, "paddingY must be > 0")
assert(cfg.paddingZ > 0, "paddingZ must be > 0")
assert(cfg.terminalPaddingFactor >= 0, "terminalPaddingFactor must be >= 0")
assert(cfg.interQuality > 0, "interQuality must be > 0")
assert(cfg.boxRefs >= 0 and cfg.boxRefs <= 10,
	"boxRefs must be an integer between 0 and 10")
ensure_dir(dirname(cfg.out))

-- Stage 1: smooth, radius-normalize, and merge all traces ---------------------
local out_base = strip_extension(cfg.out)
local merged_swc = out_base .. "_merged_outer_radius.swc"
local merged = assert(io.open(merged_swc, "w"), "Cannot create merged SWC.")
local next_id = 1
local outer_radius = cfg.radius + cfg.membraneThickness

if cfg.centerlineScale ~= 1.0 then
	print(string.format("WARNING: -centerlineScale %.6g multiplies the trace length; use 1 "
		.. "for SWC coordinates already in micrometres (1 um voxels).", cfg.centerlineScale))
end

-- Split a (merged, multi-root) SWC into its root-to-terminal chains.
local function split_chains(points)
	local by_id, children, roots = {}, {}, {}
	for _, p in ipairs(points) do
		local id, parent = tonumber(p[1]), tonumber(p[7])
		by_id[id] = p
		if parent == -1 then roots[#roots + 1] = id
		else
			assert(children[parent] == nil, "Branched SWC traces are not supported.")
			children[parent] = id
		end
	end
	local chains = {}
	for _, r in ipairs(roots) do
		local chain, id = {}, r
		while id do chain[#chain + 1] = by_id[id]; id = children[id] end
		chains[#chains + 1] = chain
	end
	return chains
end

local trace_sets = {}
local num_nephrons = #swcs
local merged_metadata = {}
if cfg.mergedSwc ~= "" then
	local chains = split_chains(read_swc(cfg.mergedSwc))
	num_nephrons = #chains
	assert(num_nephrons >= 2, "-mergedSwc must contain at least two root-to-terminal trees.")
	for i, chain in ipairs(chains) do
		trace_sets[i] = chain
	end
	local f = assert(io.open(cfg.mergedSwc, "r"))
	for line in f:lines() do
		if line:match("^%s*# UG4_") then merged_metadata[#merged_metadata + 1] = line end
	end
	f:close()
	for _, line in ipairs(merged_metadata) do merged:write(line, "\n") end
	print("Using merged SWC " .. cfg.mergedSwc .. " (" .. num_nephrons .. " trees); "
		.. "no smoothing or terminal extension.")
	cfg.extendTerminals = false
else
	for trace_index, raw_swc in ipairs(swcs) do
		local ready_swc = out_base .. "_trace_" .. trace_index .. "_mesh_ready.swc"
		ensure_dir(dirname(ready_swc))
		smooth_swc(raw_swc, ready_swc)
		trace_sets[trace_index] = read_swc(ready_swc)
	end
end
local clearance_traces = {}
for trace_index, points in ipairs(trace_sets) do
	clearance_traces[trace_index] = trace_points(points, cfg.centerlineScale)
end
audit_clearance(clearance_traces, outer_radius)

for trace_index, points in ipairs(trace_sets) do
	audit_trace_curvature(points, "trace " .. trace_index, outer_radius)
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
			"Missing parent in trace " .. trace_index)
		local point_type = new_parent == -1 and 1
			or (tonumber(point[2]) == 0 and 3 or tonumber(point[2]))
		merged:write(table.concat({
			id_map[tonumber(point[1])], point_type,
			string.format("%.17g", assert(tonumber(point[3])) * cfg.centerlineScale * MESH_LENGTH_SCALE),
			string.format("%.17g", assert(tonumber(point[4])) * cfg.centerlineScale * MESH_LENGTH_SCALE),
			string.format("%.17g", assert(tonumber(point[5])) * cfg.centerlineScale * MESH_LENGTH_SCALE),
			string.format("%.17g", outer_radius * MESH_LENGTH_SCALE), new_parent
		}, " "), "\n")
	end
end
merged:close()

local import_swc = merged_swc
if cfg.extendTerminals then
	local extended_swc = out_base .. "_outerwall_extended.swc"
	local extension_command = table.concat({
		"python3", shell_quote(cfg.terminalExtender), shell_quote(merged_swc),
		"-o", shell_quote(extended_swc),
		"--padding", tostring(cfg.padding * MESH_LENGTH_SCALE),
		"--padding-x", tostring(cfg.paddingX * MESH_LENGTH_SCALE),
		"--padding-y", tostring(cfg.paddingY * MESH_LENGTH_SCALE),
		"--padding-z", tostring(cfg.paddingZ * MESH_LENGTH_SCALE),
		"--spacing", tostring(cfg.terminalSpacing * MESH_LENGTH_SCALE),
		"--box-refs", tostring(cfg.boxRefs),
		"--min-padding-factor", tostring(cfg.terminalPaddingFactor),
		"--min-port-separation", tostring(cfg.minPortSeparation * MESH_LENGTH_SCALE),
		"--relative-tolerance", tostring(cfg.coordinateRelTol),
		"--direction-tolerance", tostring(cfg.directionTolerance)
	}, " ")
	local extension_ok = os.execute(extension_command)
	assert(extension_ok == true or extension_ok == 0,
		"Terminal-to-OuterWall SWC extension failed")
	import_swc = extended_swc
end

-- Stage 2: create all nested tubes plus one shared Inter/OuterWall domain -----
local out_ugx = out_base .. ".ugx"
local lumen_scale = cfg.radius / outer_radius
configure_adaptive_nephron_meshing(
	cfg.adaptiveAnisotropy, cfg.bendAnisotropy, cfg.bendChordErrorFactor,
	cfg.gradedInter, cfg.interNearEdgeFactor, cfg.interNearDistanceFactor,
	cfg.interFarDistanceFactor, cfg.interFarEdgeFactor,
	cfg.interMesh == "prism", cfg.prismInterThickness * MESH_LENGTH_SCALE)
configure_nephron_inter_fade(cfg.interFadeDistance * MESH_LENGTH_SCALE)
configure_nephron_geometry_policy(
	cfg.coordinateRelTol, cfg.containmentRelTol, cfg.parameterTol,
	cfg.mergeRadiusFactor, cfg.terminalSnapRadiusFactor,
	cfg.straightnessRadiusFactor, cfg.angularSpacingToleranceFactor,
	cfg.interMinLayers, cfg.interMaxLayers,
	cfg.interLayersPerSurfaceEdge, cfg.interRadialGrowth)
import_multiple_nephrons_with_membrane_and_inter_from_swc(
	import_swc, out_ugx, num_nephrons, lumen_scale, cfg.anisotropy,
	cfg.numRefs, cfg.padding * MESH_LENGTH_SCALE, cfg.interQuality, cfg.boxRefs,
	cfg.ogrid, cfg.coarseLumenCenterOnly, cfg.crossOnly)

-- Stage 3: fail early if refinement metadata was not serialized ---------------
assert(file_contains(out_ugx, "<projector type=\"neurite\""),
	"UGX is missing its NeuriteProjector.")
assert(file_contains(out_ugx, "<vertex_attachment name=\"npSurfParams\""),
	"UGX is missing npSurfParams.")
validate_ugx(out_ugx)
for level = 1, cfg.numRefs do
	validate_ugx(out_base .. "_refined_" .. level .. ".ugx")
end

-- Stage 4: concise result summary --------------------------------------------
local vertices, edges, faces, volumes = count_mesh_elements(out_ugx)
print("Nephrons: " .. num_nephrons)
print("Vertices: " .. vertices)
print("Edges: " .. edges)
print("Faces: " .. faces)
print("Volumes: " .. volumes)
print("Lumen diameter (um): " .. (2.0 * cfg.radius))
print("Membrane thickness (um): " .. cfg.membraneThickness)
print("Outer radius (um): " .. outer_radius)
print("Adaptive anisotropy straight/bend: " .. cfg.anisotropy .. " / "
	.. cfg.bendAnisotropy)
print("Graded Inter: " .. tostring(cfg.gradedInter) .. ", edge factors "
	.. cfg.interNearEdgeFactor .. " -> " .. cfg.interFarEdgeFactor
	.. ", distance factors " .. cfg.interNearDistanceFactor .. " -> "
	.. cfg.interFarDistanceFactor)
print("Inter mesh topology: " .. cfg.interMesh)
print("Structured near-Inter thickness (um): " .. cfg.prismInterThickness)
print("Centerline scale: " .. cfg.centerlineScale)
print("SWC resample spacing (um): " .. cfg.resampleSpacing)
print("Terminal extension spacing (um): " .. cfg.terminalSpacing)
print("OuterWall padding X/Y/Z (um): " .. cfg.paddingX .. " / "
	.. cfg.paddingY .. " / " .. cfg.paddingZ)
print("Mesh coordinate unit: mm")
print("Terminal-to-OuterWall extension: " .. tostring(cfg.extendTerminals))
print("Output: " .. out_ugx)
