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

-- Build a projector-backed neurite surface directly from the specified SWC tree.
--
-- This is a copy of swc_to_neurite_projector.lua tailored to the
-- Nephron traces dataset and a fixed radius of 1.0.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

-- the script directory is used to resolve relative paths for the default SWC and output UGX
-- the input paramaters with call out to the terminal to override the defaults
local script_dir = GENERATION_DIR

local cfg = {}
-- SWC coordinates and CLI lengths are in um; UGX/projector coordinates are mm.
local MESH_LENGTH_SCALE = 1e-3
cfg.swc = util.GetParam("-swc", MODEL_ROOT .. "/Nephron traces/3rd_trace_243.swc")
local input_name = cfg.swc:match("([^/]+)$") or cfg.swc
local input_stem = input_name:gsub("%.[^%.]+$", "")
cfg.out = util.GetParam("-out", MODEL_ROOT .. "/runs/" .. input_stem .. ".ugx")
cfg.smoother = util.GetParam("-smoother",
	script_dir .. "/quadratic_interpolate_swc.py")
cfg.readySwc = util.GetParam("-readySwc", "")
cfg.extendTerminals = util.GetParamNumber("-extendTerminals", 1) ~= 0
cfg.terminalExtender = util.GetParam("-terminalExtender",
	script_dir .. "/extend_swc_terminals_to_outerwall.py")
cfg.terminalPaddingFactor = util.GetParamNumber("-terminalPaddingFactor", 25.0)
-- Scale the centerline independently from its tube radii.  This is required
-- when a requested outer radius is too large for the curvature/clearance of
-- the original trace.  A value of 1 preserves the measured centerline.
cfg.centerlineScale = util.GetParamNumber("-centerlineScale", 1.0)
cfg.smoothOrder = math.floor(util.GetParamNumber("-smoothOrder", 3))
cfg.samplesPerCorner = math.floor(util.GetParamNumber("-samplesPerCorner", 30))
-- Cubic interpolation can reintroduce 20--25 degree turns even after the raw
-- controls were softened at 25 degrees.  A 10-degree control threshold keeps
-- the coarse axial sweep valid when refinement uses the exact projector.
cfg.preSmoothAngle = util.GetParamNumber("-preSmoothAngle", 10.0)
cfg.preSmoothStrength = util.GetParamNumber("-preSmoothStrength", 0.85)
cfg.preSmoothPasses = math.floor(util.GetParamNumber("-preSmoothPasses", 50))
cfg.radius = util.GetParamNumber("-radius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
local requested_outer_radius = cfg.radius + cfg.membraneThickness
-- Smoothing operates on the original SWC coordinates.  Keep its sampling tied
-- to that source trace; centerlineScale is applied afterwards.
cfg.resampleSpacing = util.GetParamNumber(
	"-resampleSpacing", 0.5)
-- Terminal-extension samples may be coarser for a large tube.  This avoids
-- thousands of nearly redundant spline sections while retaining the verified
-- 0.5 um spacing for the 1.0/0.2 um reference geometry.
cfg.terminalSpacing = util.GetParamNumber("-terminalSpacing",
	math.max(cfg.resampleSpacing * cfg.centerlineScale,
		0.25 * requested_outer_radius))
cfg.skipGeometryPreflight = util.GetParamNumber("-skipGeometryPreflight", 0) ~= 0
-- Keep long axial cells only where the spline is nearly straight. Around
-- bends, recursively add stations until the axial/circumferential ratio is
-- bendAnisotropy (1 by default).
cfg.anisotropy = util.GetParamNumber("-anisotropy", 4.0)
cfg.adaptiveAnisotropy = util.GetParamNumber("-adaptiveAnisotropy", 1) ~= 0
cfg.bendAnisotropy = util.GetParamNumber("-bendAnisotropy", 1.0)
cfg.bendChordErrorFactor = util.GetParamNumber("-bendChordErrorFactor", 0.02)
-- Grade Inter from the Basolateral surface into the outer box. Every factor
-- is relative to the Basolateral boundary-edge scale, not a fixed length.
cfg.gradedInter = util.GetParamNumber("-gradedInter", 1) ~= 0
cfg.interNearEdgeFactor = util.GetParamNumber("-interNearEdgeFactor", 2.0)
cfg.interNearDistanceFactor = util.GetParamNumber("-interNearDistanceFactor", 1.5)
cfg.interFarDistanceFactor = util.GetParamNumber("-interFarDistanceFactor", 6.0)
cfg.interFarEdgeFactor = util.GetParamNumber("-interFarEdgeFactor", 6.0)
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
-- The non-terminal box clearance is tied both to the centerline and tube.
cfg.padding = util.GetParamNumber("-padding",
	math.max(10.0 * cfg.centerlineScale, 2.0 * requested_outer_radius))
-- Axis-specific clearances inherit the common padding unless overridden.
-- They control the fixed box metadata created by the terminal extender.
cfg.paddingX = util.GetParamNumber("-paddingX", cfg.padding)
cfg.paddingY = util.GetParamNumber("-paddingY", cfg.padding)
cfg.paddingZ = util.GetParamNumber("-paddingZ", cfg.padding)
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

local function write_radius_override_swc(src, dst, radius, centerline_scale)
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
				string.format("%.17g", assert(tonumber(tokens[3])) * centerline_scale * MESH_LENGTH_SCALE),
				string.format("%.17g", assert(tonumber(tokens[4])) * centerline_scale * MESH_LENGTH_SCALE),
				string.format("%.17g", assert(tonumber(tokens[5])) * centerline_scale * MESH_LENGTH_SCALE),
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

local function geometry_preflight(path, outer_radius_um)
	local points = {}
	local file = assert(io.open(path, "r"), "Cannot audit SWC: " .. path)
	for line in file:lines() do
		local payload = line:gsub("#.*$", "")
		local fields = {}
		for token in payload:gmatch("%S+") do fields[#fields + 1] = token end
		if #fields > 0 then
			assert(#fields == 7, "SWC line must have 7 columns: " .. line)
			points[#points + 1] = {
				x = assert(tonumber(fields[3])) * cfg.centerlineScale,
				y = assert(tonumber(fields[4])) * cfg.centerlineScale,
				z = assert(tonumber(fields[5])) * cfg.centerlineScale
			}
		end
	end
	file:close()
	assert(#points >= 3, "Geometry preflight needs at least three SWC points")

	local function distance(a, b)
		local dx, dy, dz = a.x-b.x, a.y-b.y, a.z-b.z
		return math.sqrt(dx*dx + dy*dy + dz*dz)
	end
	local cumulative = {0.0}
	for i = 2, #points do
		cumulative[i] = cumulative[i-1] + distance(points[i-1], points[i])
	end

	local min_curvature_radius = math.huge
	for i = 2, #points-1 do
		local a, b, c = points[i-1], points[i], points[i+1]
		local ab = {x=b.x-a.x, y=b.y-a.y, z=b.z-a.z}
		local ac = {x=c.x-a.x, y=c.y-a.y, z=c.z-a.z}
		local cross = {
			x=ab.y*ac.z-ab.z*ac.y,
			y=ab.z*ac.x-ab.x*ac.z,
			z=ab.x*ac.y-ab.y*ac.x
		}
		local twice_area = math.sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z)
		if twice_area > 1e-12 * cfg.centerlineScale * cfg.centerlineScale then
			local circumradius = distance(a,b) * distance(b,c) * distance(a,c)
				/ (2.0 * twice_area)
			min_curvature_radius = math.min(min_curvature_radius, circumradius)
		end
	end

	local outer_radius = outer_radius_um
	local exclusion_arclength = math.max(4.0 * outer_radius,
		4.0 * cfg.resampleSpacing * cfg.centerlineScale)
	local min_nonlocal_distance = math.huge
	for i = 1, #points do
		for j = i+1, #points do
			if cumulative[j] - cumulative[i] > exclusion_arclength then
				min_nonlocal_distance = math.min(min_nonlocal_distance,
					distance(points[i], points[j]))
			end
		end
	end

	print(string.format(
		"Geometry preflight (um): outerRadius=%.9g, minCurvatureRadius=%.9g, minNonlocalCenterlineDistance=%.9g",
		outer_radius, min_curvature_radius, min_nonlocal_distance))
	if cfg.skipGeometryPreflight then return end
	assert(outer_radius < 0.9 * min_curvature_radius,
		string.format("Requested outer radius %.6g um is incompatible with centerline curvature (minimum radius %.6g um). Reduce radius/thickness or use -centerlineScale > %.6g.",
			outer_radius, min_curvature_radius,
			outer_radius / (0.9 * min_curvature_radius) * cfg.centerlineScale))
	assert(2.0 * outer_radius < 0.9 * min_nonlocal_distance,
		string.format("Requested outer diameter %.6g um exceeds non-local centerline clearance %.6g um. Reduce radius/thickness or enlarge the centerline.",
			2.0 * outer_radius, min_nonlocal_distance))
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
assert(cfg.anisotropy >= 1.0, "anisotropy must be >= 1")
assert(cfg.bendAnisotropy >= 1.0 and cfg.bendAnisotropy <= cfg.anisotropy,
	"bendAnisotropy must be between 1 and anisotropy")
assert_positive("bendChordErrorFactor", cfg.bendChordErrorFactor)
assert_positive("interNearEdgeFactor", cfg.interNearEdgeFactor)
assert(cfg.interFarEdgeFactor >= cfg.interNearEdgeFactor,
	"interFarEdgeFactor must be >= interNearEdgeFactor")
assert(cfg.interNearDistanceFactor >= 0.0
	and cfg.interFarDistanceFactor > cfg.interNearDistanceFactor,
	"Inter distance factors require 0 <= near < far")
assert(cfg.ogrid >= 4, "ogrid must be an integer >= 4")
assert_positive("padding", cfg.padding)
assert_positive("paddingX", cfg.paddingX)
assert_positive("paddingY", cfg.paddingY)
assert_positive("paddingZ", cfg.paddingZ)
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
assert_positive("centerlineScale", cfg.centerlineScale)
assert_positive("terminalSpacing", cfg.terminalSpacing)
assert_positive("membraneThickness", cfg.membraneThickness)
assert(cfg.numRefs >= 0, "numRefs must be >= 0")
if cfg.radius == 0.0 then
	error("radius must be > 0, or < 0 to keep the SWC radii")
end

ensure_dir(dirname(cfg.out))

local ready_swc = cfg.readySwc
if ready_swc == "" then
	-- Output-specific preprocessing prevents one parameter sweep from silently
	-- overwriting the ready SWC used by another build.
	ready_swc = strip_extension(cfg.out) .. "_mesh_ready.swc"
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
	geometry_preflight(ready_swc, outer_radius)
	write_radius_override_swc(ready_swc, import_swc, outer_radius,
		cfg.centerlineScale)
else
	error("The membrane test currently requires -radius > 0 so it can apply a constant physical membrane thickness.")
end
local lumen_scale = cfg.radius / outer_radius

configure_adaptive_nephron_meshing(
	cfg.adaptiveAnisotropy, cfg.bendAnisotropy, cfg.bendChordErrorFactor,
	cfg.gradedInter, cfg.interNearEdgeFactor, cfg.interNearDistanceFactor,
	cfg.interFarDistanceFactor, cfg.interFarEdgeFactor)

if cfg.extendTerminals then
	local extended_swc = strip_extension(cfg.out) .. "_outerwall_extended.swc"
	local extension_command = table.concat({
		"python3", shell_quote(cfg.terminalExtender), shell_quote(import_swc),
		"-o", shell_quote(extended_swc),
		"--padding", tostring(cfg.padding * MESH_LENGTH_SCALE),
		"--padding-x", tostring(cfg.paddingX * MESH_LENGTH_SCALE),
		"--padding-y", tostring(cfg.paddingY * MESH_LENGTH_SCALE),
		"--padding-z", tostring(cfg.paddingZ * MESH_LENGTH_SCALE),
		"--spacing", tostring(cfg.terminalSpacing * MESH_LENGTH_SCALE),
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
print("Membrane thickness (um): " .. cfg.membraneThickness)
print("Outer radius (um): " .. outer_radius)
print("Adaptive anisotropy straight/bend: " .. cfg.anisotropy .. " / "
	.. cfg.bendAnisotropy)
print("Graded Inter: " .. tostring(cfg.gradedInter) .. ", edge factors "
	.. cfg.interNearEdgeFactor .. " -> " .. cfg.interFarEdgeFactor
	.. ", distance factors " .. cfg.interNearDistanceFactor .. " -> "
	.. cfg.interFarDistanceFactor)
print("Centerline scale: " .. cfg.centerlineScale)
print("SWC resample spacing (um): " .. cfg.resampleSpacing)
print("Terminal extension spacing (um): " .. cfg.terminalSpacing)
print("OuterWall padding X/Y/Z (um): " .. cfg.paddingX .. " / "
	.. cfg.paddingY .. " / " .. cfg.paddingZ)
print("Mesh coordinate unit: mm")
print("Terminal-to-OuterWall extension: " .. tostring(cfg.extendTerminals))
print("Output: " .. out_ugx)
