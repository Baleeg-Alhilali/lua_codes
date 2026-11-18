--------------------------------------------------------------------------------
-- Prepare Long_tube2 for the Stokes--Darcy projector control
--------------------------------------------------------------------------------
-- Some Long_tube2 exports retain an internal lumen cross-section in Inlet,
-- with or without also marking the physical z=250 cap. This script removes
-- every old Inlet face, then assigns only the positive-z lumen cap to Inlet.
-- It changes face subset assignments only and preserves the source mesh.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh"})

local source = "runs/Long_tube2.ugx"
local target = "runs/Long_tube2_flow_control.ugx"
local eps = 1.0e-6

mesh = Mesh()
assert(LoadMesh(mesh, source), "Could not load " .. source)
local sh = mesh:subset_handler()
local lumenSubset = sh:get_subset_index("Lumen")
local inletSubset = sh:get_subset_index("Inlet")
assert(lumenSubset >= 0, "The source mesh has no Lumen subset")
assert(inletSubset >= 0, "The source mesh has no Inlet subset")

-- Remove the invalid internal inlet label.
DeselectAllFaces(mesh)
SelectSubset(mesh, inletSubset, false, false, true, false)
AssignSubset(mesh, lumenSubset, false, false, true, false)

-- Assign the physical positive-z lumen cap as the inlet.
DeselectAllFaces(mesh)
SelectFacesInCylinder(
    mesh,
    Vec3d(0.0, 0.0, 250.0 - eps),
    Vec3d(0.0, 0.0, 250.0 + eps),
    10.001
)
RestrictSelectionToSubset(mesh, lumenSubset)
AssignSubset(mesh, inletSubset, false, false, true, false)

DeselectAllFaces(mesh)
assert(SaveMesh(mesh, target), "Could not write " .. target)
print("Wrote boundary-corrected flow-control mesh: " .. target)
