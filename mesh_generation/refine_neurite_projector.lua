-- Resolve default grids from Model while running inside mesh_generation.
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

-- Refine an existing neurite-projector UGX while preserving projector metadata.
--
-- Input must contain both:
--   * <projector type="neurite">
--   * <vertex_attachment name="npSurfParams">
--
-- If npSurfParams is missing, the neurite projector cannot recover the
-- per-vertex axial/angular/radial coordinates needed to project new vertices.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh", "neuro_collection"})
InitUG(3, AlgebraType("CPU", 1))

local cfg = {}
cfg.grid = util.GetParam("-grid", MODEL_ROOT .. "/runs/sample_tree_projector.ugx")
cfg.out = util.GetParam("-out", MODEL_ROOT .. "/runs/sample_tree_projector_refined.ugx")
cfg.refs = math.floor(util.GetParamNumber("-refs", 1))
cfg.validate = util.GetParamNumber("-validate", 1) ~= 0
cfg.validator = util.GetParam("-validator", GENERATION_DIR .. "/validate_ugx_mesh.py")

local function shell_quote(s)
	return "'" .. tostring(s):gsub("'", "'\\''") .. "'"
end

local function dirname(path)
	local dir = path:match("^(.*)/[^/]*$")
	if dir == nil or dir == "" then
		return "."
	end
	return dir
end

local function ensure_dir(path)
	if path ~= "." then
		local ok = os.execute("mkdir -p " .. shell_quote(path))
		assert(ok == true or ok == 0, "Could not create directory: " .. path)
	end
end

local function file_contains(path, needle)
	local file = assert(io.open(path, "r"), "Cannot open file: " .. path)
	local text = file:read("*a")
	file:close()
	return text:find(needle, 1, true) ~= nil
end

local function assert_neurite_projector_ugx(path)
	assert(file_contains(path, "<projector type=\"neurite\""),
		"UGX is missing a serialized NeuriteProjector: " .. path)
	assert(file_contains(path, "<vertex_attachment name=\"npSurfParams\""),
		"UGX is missing npSurfParams. Regenerate it from the SWC file, or use an earlier UGX that still contains npSurfParams: " .. path)
end

assert(cfg.grid ~= "", "Pass an input UGX with -grid <path>")
assert(cfg.out ~= "", "Pass an output path with -out <path>")
assert(cfg.refs >= 0, "refs must be >= 0")

ensure_dir(dirname(cfg.out))
assert_neurite_projector_ugx(cfg.grid)

print("Neurite projector refine")
print("  input UGX = " .. cfg.grid)
print("  output    = " .. cfg.out)
print("  refs      = " .. cfg.refs)

local mesh = Mesh()
assert(LoadMesh(mesh, cfg.grid), "LoadMesh failed: " .. cfg.grid)

for i = 1, cfg.refs do
	SelectAll(mesh)
	print("  projected refine " .. i .. " / " .. cfg.refs)
	Refine(mesh)

	-- Projected refinement can leave detached child vertices at constrained
	-- terminal interfaces.  They are not part of any edge, face, or volume and
	-- must not be serialized or allowed to become parents of the next level.
	ClearSelection(mesh)
	local orphanVertices = SelectUnconnectedVertices(mesh, true, true, true)
	if orphanVertices > 0 then
		EraseSelectedElements(mesh, true, true, true)
	end

	-- Projection changes child coordinates after their connectivity has been
	-- created.  A non-degenerate tetrahedron may consequently have the opposite
	-- vertex ordering.  Normalize every volume after every level, exactly as the
	-- Neuro Collection generator does before writing its flat level files.
	ClearSelection(mesh)
	SelectAllVolumes(mesh)
	local corrected = FixVolumeOrientation(mesh)

	-- Check orientation determinacy after repair.  A cell that has the same
	-- orientation before and after a flip is collapsed/pathological and is a
	-- hard mesh-generation failure, never something to pass to a solver.
	ClearSelection(mesh)
	local unorientable = SelectUnorientableVolumes(mesh)
	assert(unorientable == 0,
		"Projected refinement level " .. i .. " created " .. unorientable
		.. " unorientable/degenerate volumes")
	ClearSelection(mesh)

	print("    removed detached vertices: " .. orphanVertices)
	print("    corrected volume orientations: " .. corrected)
	print("    unorientable/degenerate volumes: " .. unorientable)
end

assert(SaveMesh(mesh, cfg.out), "SaveMesh failed: " .. cfg.out)
assert_neurite_projector_ugx(cfg.out)
if cfg.validate then
	print("Validating refined mesh: " .. cfg.out)
	local ok = os.execute("python3 " .. shell_quote(cfg.validator)
		.. " " .. shell_quote(cfg.out))
	assert(ok == true or ok == 0,
		"Mesh validation failed; refusing refined output: " .. cfg.out)
end

print("Done. Output still contains NeuriteProjector and npSurfParams.")
