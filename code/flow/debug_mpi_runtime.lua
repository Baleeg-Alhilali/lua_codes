ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
InitUG(3, AlgebraType("CPU", 1))
dom = Domain()
LoadDomain(dom, "runs/straight_nephron_no_bend_plain.ugx")
print("pre-distribution multigrid consistency: " .. tostring(CheckMultiGridConsistency(dom:grid())))
print("pre-distribution hanging-node consistency: " .. tostring(CheckHangingNodeConsistency(dom:grid())))
balancer.firstDistLvl = 0
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, 0)
print("after redistribution rank " .. ProcRank())
print("post-distribution multigrid consistency: " .. tostring(CheckMultiGridConsistency(dom:grid())))
print("post-distribution hanging-node consistency: " .. tostring(CheckHangingNodeConsistency(dom:grid())))
if true then
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("p_darcy_scaled", "Lagrange", 1,
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet,TripleJunction")
approxSpace:init_levels()
approxSpace:init_top_surface()
print("MPI mesh rank " .. ProcRank() .. " of " .. NumProcs())
end
