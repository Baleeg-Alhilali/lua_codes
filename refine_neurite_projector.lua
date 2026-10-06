-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

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
cfg.grid = util.GetParam("-grid", "runs/sample_tree_projector.ugx")
cfg.out = util.GetParam("-out", "runs/sample_tree_projector_refined.ugx")
cfg.refs = math.floor(util.GetParamNumber("-refs", 1))

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
end

assert(SaveMesh(mesh, cfg.out), "SaveMesh failed: " .. cfg.out)
assert_neurite_projector_ugx(cfg.out)

print("Done. Output still contains NeuriteProjector and npSurfParams.")
