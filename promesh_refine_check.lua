-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

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

local grid = util.GetParam("-grid", "runs/sample_tree_projector.ugx")
local out = util.GetParam("-out", "runs/promesh_refined.ugx")

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
