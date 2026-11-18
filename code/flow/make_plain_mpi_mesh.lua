ug_load_script("ug_util.lua")
InitUG(3, AlgebraType("CPU", 1))

local dom = Domain()
LoadDomain(dom, "runs/straight_nephron_no_bend.ugx")
print("original multigrid consistency: " .. tostring(CheckMultiGridConsistency(dom:grid())))
print("original hanging-node consistency: " .. tostring(CheckHangingNodeConsistency(dom:grid())))
dom:set_refinement_projector(
    ProjectionHandler(dom:geometry3d(), dom:subset_handler()))
SaveDomain(dom, "runs/straight_nephron_no_bend_plain.ugx")
print("Saved projector-free MPI diagnostic mesh.")
