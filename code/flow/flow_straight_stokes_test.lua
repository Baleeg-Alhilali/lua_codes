--------------------------------------------------------------------------------
-- Diagnostic summary
--------------------------------------------------------------------------------
-- Isolates the lumen Stokes FV1 operator and its pressure stabilization from
-- Darcy and PressureTraceConstraint. Apical is temporarily a no-slip wall.
-- A direct solve is used so a refinement-born pressure mode cannot be blamed
-- on GMG. Compare level 0 and level 1 lumen pressure outputs.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------

local scriptPath = debug.getinfo(1, "S").source
if scriptPath:sub(1, 1) == "@" then scriptPath = scriptPath:sub(2) end
if scriptPath:sub(1, 1) ~= "/" then
    scriptPath = CurrentWorkingDirectory() .. "/" .. scriptPath
end
local MODEL_ROOT = assert(scriptPath:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. scriptPath)
ChangeDirectory(MODEL_ROOT)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
AssertPluginsLoaded({"NavierStokes", "neuro_collection", "SuperLU6"})

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3
numRefs = util.GetParamNumber("-numRefs", 0, "number of uniform refinements")
gridName = "runs/straight_nephron_no_bend.ugx"
outputPrefix = util.GetParam(
    "-output",
    "Results_final/debug_old_mesh_stokes_fields_cor_ref" .. numRefs,
    "VTK output prefix"
)

vmax = 3.0                         -- [mm/s]
viscosity = 1.0                    -- [mm^2/s]
inletRadius = 1.5e-2               -- [mm]
inletCenterX = 0.0                 -- [mm]
inletCenterY = 0.0                 -- [mm]

assert(numRefs >= 0, "numRefs must be non-negative")
InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"}),
    "Stokes control mesh is missing a required subset")

balancer.firstDistLvl = numRefs
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, numRefs)
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
local stokesSubsets = "Lumen,Apical,Inlet,Outlet,TripleJunction"
approxSpace:add_fct({"u", "v", "w", "p_stokes"},
    "Lagrange", 1, {"Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"})
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Stokes discretization and boundary conditions
--------------------------------------------------------------------------------

function inletVelocity(x, y, z, t)
    local dx, dy = x - inletCenterX, y - inletCenterY
    local r2 = dx * dx + dy * dy
    local profile = math.max(0.0, 1.0 - r2 / (inletRadius * inletRadius))
    if r2 >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return 0.0, 0.0, vmax * profile
end

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
stokesEq:set_stokes(true)
stokesEq:set_laplace(false)
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_upwind("lps")
stokesEq:set_peclet_blend(false)
stokesEq:set_stabilization("fields", "cor")
stokesEq:set_pac_upwind(false)

inlet = NavierStokesInflow(stokesEq)
inlet:add("inletVelocity", "Inlet")
outlet = NavierStokesNoNormalStressOutflow(stokesEq)
outlet:add("Outlet,Apical")
wall = NavierStokesWall(stokesEq)
wall:add("TripleJunction")
bjs = NavierStokesWSBCFV1(stokesEq)
bjs:add("Apical")
bjs:set_sliding_factor(2.0e2)
bjs:set_sliding_limit(0.0)
bjs:set_normal_factor(1.0e3)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(stokesEq)
domainDisc:add(inlet)
domainDisc:add(outlet)
domainDisc:add(wall)
domainDisc:add(bjs)

--------------------------------------------------------------------------------
-- Uzawa-GMG diagnostic solver
--------------------------------------------------------------------------------

matrix = AssembledLinearOperator(domainDisc)
rhs = GridFunction(approxSpace)
solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)
domainDisc:assemble_linear(matrix, rhs)

mainInverse = ILU()
mainInverse:set_damp(1.0)
mainInverse:set_inversion_eps(1.0e-30)
pressureInverse = ILU()
pressureInverse:set_damp(1.0)
pressureInverse:set_inversion_eps(1.0e-30)
smoother = UzawaBase({"p_stokes"})
smoother:set_forward_iter(mainInverse)
smoother:set_schur_iter(pressureInverse)
smoother:set_backward_iter(mainInverse)

schurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
schurMassEq:set_mass_scale(1.0)
schurMassEq:set_diffusion(0.0)
schurMassDisc = DomainDiscretization(approxSpace)
schurMassDisc:add(schurMassEq)
schurMassOp = AssembledLinearOperator(schurMassDisc)
-- No mass update: use the stabilized FV1 pressure block as assembled.

fullSmoother = ILU()
fullSmoother:set_damp(1.0)
fullSmoother:set_inversion_eps(1.0e-30)
fullSmoother:enable_consistent_interfaces(true)

solver = util.solver.CreateLinearSolver({
    type = "bicgstab",
    precond = {
        type = "gmg", approxSpace = approxSpace,
        discretization = domainDisc, baseLevel = 0,
        baseSolver = "superlu", smoother = fullSmoother,
        cycle = "V", preSmooth = 2, postSmooth = 2, rap = false,
        gatheredBaseSolverIfAmbiguous = true,
        transfer = {type = "std", restrictionDamp = 1.0,
            prolongationDamp = 1.0, enableP1LagrangeOptimization = true}
    },
    convCheck = {type = "standard", iterations = 200,
        absolute = 1.0e-12, reduction = 1.0e-8, verbose = true}
})
solver:init(matrix, solution)
assert(solver:apply(solution, rhs), "Stokes control GMG solve failed")
solution:set_consistent_storage_type()

defect = rhs:clone()
VecAssign(defect, rhs)
matrix:apply_sub(defect, solution)
print(string.format("Stokes control true defect: %.12e", VecNorm(defect)))

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

output = VTKOutput()
output:clear_selection()
output:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
output:select_nodal("p_stokes", "stokes_pressure_kinematic")
output:print_subsets(outputPrefix .. "_lumen", solution, "Lumen", 0, 0.0)

print("Stokes refinement control output prefix: " .. outputPrefix)
