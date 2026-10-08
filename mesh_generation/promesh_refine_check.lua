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

-- Reproduce ProMesh's LoadMesh + Refine path on a neurite-projector UGX.
--
-- This is intentionally close to the GUI/tool path: ProMesh LoadMesh reads
-- geometry into a Mesh, then Refine uses Mesh:projection_handler().

ug_load_script("ug_util.lua")

local loadNeuro = util.GetParamNumber("-loadNeuro", 1)
if loadNeuro == 0 then
	AssertPluginsLoaded({"ProMesh"})
else
	AssertPluginsLoaded({"ProMesh", "neuro_collection"})
end
InitUG(3, AlgebraType("CPU", 1))

local grid = util.GetParam("-grid", MODEL_ROOT .. "/runs/sample_tree_projector.ugx")
local out = util.GetParam("-out", MODEL_ROOT .. "/runs/promesh_refined.ugx")

print("ProMesh refine check")
print("  grid = " .. grid)
print("  out  = " .. out)
print("  neuro_collection = " .. (loadNeuro == 0 and "not loaded by script" or "loaded"))

local mesh = Mesh()
assert(LoadMesh(mesh, grid), "LoadMesh failed: " .. grid)
SelectAll(mesh)
Refine(mesh)
SaveMesh(mesh, out)

print("Done.")
