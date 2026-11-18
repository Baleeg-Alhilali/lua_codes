--------------------------------------------------------------------------------
-- Extract a regular Lumen-only companion for high-refinement flow controls
--------------------------------------------------------------------------------
-- The source mesh is preserved.  Removing the Membrane and Inter volumes
-- before refinement avoids paying their much larger refinement cost in a
-- diagnostic whose approximation space and equations live only on Lumen.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh"})

local source = util.GetParam("-grid", "runs/Long_tube_V3.ugx",
    "whole-domain source UGX")
local target = util.GetParam("-out", "runs/Long_tube_V3_lumen.ugx",
    "new Lumen-only companion UGX")

mesh = Mesh()
assert(LoadMesh(mesh, source), "Could not load source mesh: " .. source)
local sh = mesh:subset_handler()
local lumenSubset = sh:get_subset_index("Lumen")
local membraneSubset = sh:get_subset_index("Membrane")
local interSubset = sh:get_subset_index("Inter")
assert(lumenSubset >= 0, "Source mesh has no Lumen subset")
assert(membraneSubset >= 0, "Source mesh has no Membrane subset")
assert(interSubset >= 0, "Source mesh has no Inter subset")

-- Select only the two unwanted volume sets.  EraseSelectedElements then
-- removes lower-dimensional objects only when they are no longer used by a
-- retained Lumen volume; the Apical, Inlet, Outlet, and TripleJunction
-- boundary objects therefore remain attached to Lumen.
DeselectAllVertices(mesh)
DeselectAllEdges(mesh)
DeselectAllFaces(mesh)
DeselectAllVolumes(mesh)
SelectSubset(mesh, membraneSubset, false, false, false, true)
SelectSubset(mesh, interSubset, false, false, false, true)
EraseSelectedElements(mesh, true, true, true)
DeselectAllVertices(mesh)
DeselectAllEdges(mesh)
DeselectAllFaces(mesh)
DeselectAllVolumes(mesh)

assert(SaveMesh(mesh, target), "Could not write Lumen companion: " .. target)
print("Wrote Lumen-only flow companion: " .. target)
