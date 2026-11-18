--------------------------------------------------------------------------------
-- Diagnostic summary
--------------------------------------------------------------------------------
-- Isolates the lumen Stokes FV1 operator and its pressure stabilization from
-- Darcy and PressureTraceConstraint. Apical is temporarily a no-slip wall.
-- The regular Lumen-only companion removes unrelated Membrane/Inter refinement
-- cost without changing the Lumen mesh, boundary subsets, or projectors.

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
gridName = "runs/Long_tube_V3_lumen.ugx"
outputPrefix = util.GetParam(
    "-output",
    "Results_final/debug_old_mesh_stokes_fields_cor_ref" .. numRefs,
    "VTK output prefix"
)

vmax = 3.0                         -- [mm/s]
viscosity = 1.0                    -- [mm^2/s]
inletRadius = 1.0e-2               -- [mm]
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
ScaleDomain(dom, 1.0e-3, 1.0e-3, 1.0e-3)
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
stokesEq:set_laplace(true)
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_upwind("lps")
stokesEq:set_peclet_blend(false)
stokesEq:set_stabilization("fields", "cor")
stokesEq:set_pac_upwind(false)

inlet = NavierStokesInflow(stokesEq)
inlet:add("inletVelocity", "Inlet")
outlet = NavierStokesNoNormalStressOutflow(stokesEq)
outlet:add("Outlet")
wall = NavierStokesWall(stokesEq)
wall:add("TripleJunction,Apical")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(stokesEq)
domainDisc:add(inlet)
domainDisc:add(outlet)
domainDisc:add(wall)

--------------------------------------------------------------------------------
-- ILU-smoothed GMG diagnostic solver
--------------------------------------------------------------------------------

matrix = AssembledLinearOperator(domainDisc)
rhs = GridFunction(approxSpace)
solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)
domainDisc:assemble_linear(matrix, rhs)

fullSmoother = ILU()
fullSmoother:set_damp(0.7)
fullSmoother:set_inversion_eps(1.0e-30)
fullSmoother:enable_consistent_interfaces(true)

solver = util.solver.CreateLinearSolver({
    type = "bicgstab",
    precond = {
        type = "gmg", approxSpace = approxSpace,
        discretization = domainDisc, baseLevel = 0,
        baseSolver = "superlu", smoother = fullSmoother,
        cycle = "V", preSmooth = 4, postSmooth = 4, rap = false,
        gatheredBaseSolverIfAmbiguous = true,
        transfer = {type = "std", restrictionDamp = 1.0,
            prolongationDamp = 1.0, enableP1LagrangeOptimization = true}
    },
    convCheck = {type = "standard", iterations = 200,
        absolute = 1.0e-12, reduction = 1.0e-12, verbose = true}
})
solver:init(matrix, solution)
assert(solver:apply(solution, rhs), "Stokes control GMG solve failed")
solution:set_consistent_storage_type()

defect = rhs:clone()
VecAssign(defect, rhs)
matrix:apply_sub(defect, solution)
trueDefectNorm = VecNorm(defect)
print(string.format("Stokes control true defect: %.12e", trueDefectNorm))
assert(trueDefectNorm <= 1.0e-11,
    "Stokes control recomputed defect exceeds 1e-11")

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

output = VTKOutput()
output:clear_selection()
output:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
output:select_nodal("p_stokes", "stokes_pressure_kinematic")
output:print_subsets(outputPrefix .. "_lumen", solution, "Lumen", 0, 0.0)

velocityData = GridFunctionVectorData(solution, "u,v,w")
inletPressureMean = Integral(solution, "p_stokes", "Inlet")
    / Integral(1.0, solution, "Inlet")
outletPressureMean = Integral(solution, "p_stokes", "Outlet")
    / Integral(1.0, solution, "Outlet")
pressureDrop = inletPressureMean - outletPressureMean
prescribedInletFlowFV1 = IntegralNormalComponentOnManifold(
    "inletVelocity", solution, "Inlet", "Lumen", 0.0, 1)
inletFlowFV1 = IntegralNormalComponentOnManifold(
    velocityData, solution, "Inlet", "Lumen", 0.0, 1)
outletFlowFV1 = IntegralNormalComponentOnManifold(
    velocityData, solution, "Outlet", "Lumen", 0.0, 1)
apicalFlowFV1 = IntegralNormalComponentOnManifold(
    velocityData, solution, "Apical", "Lumen", 0.0, 1)
flowBalanceFV1 = prescribedInletFlowFV1 + outletFlowFV1 + apicalFlowFV1
relativeFlowBalanceFV1 = math.abs(flowBalanceFV1)
    / math.max(math.abs(prescribedInletFlowFV1), 1.0e-30)

inletFlowHighOrder = IntegralNormalComponentOnManifold(
    velocityData, solution, "Inlet", "Lumen", 0.0, 5)
outletFlowHighOrder = IntegralNormalComponentOnManifold(
    velocityData, solution, "Outlet", "Lumen", 0.0, 5)
apicalFlowHighOrder = IntegralNormalComponentOnManifold(
    velocityData, solution, "Apical", "Lumen", 0.0, 5)
flowBalanceHighOrder = inletFlowHighOrder + outletFlowHighOrder
    + apicalFlowHighOrder
relativeFlowBalanceHighOrder = math.abs(flowBalanceHighOrder)
    / math.max(math.abs(inletFlowHighOrder), 1.0e-30)
fluxBalanceTolerance = 1.0e-2
exactCircularInletFlow = -0.5 * math.pi * inletRadius * inletRadius * vmax
inletIntegrationError = math.abs(
    (prescribedInletFlowFV1 - exactCircularInletFlow) / exactCircularInletFlow)
print(string.format("Mean inlet pressure:         %.12e mm^2/s^2", inletPressureMean))
print(string.format("Mean outlet pressure:        %.12e mm^2/s^2", outletPressureMean))
print(string.format("Mean pressure drop:          %.12e mm^2/s^2", pressureDrop))
print("Flux convention: positive means outward from Lumen.")
print(string.format("Exact circular inlet flux:   %.12e mm^3/s", exactCircularInletFlow))
print(string.format("FV1 prescribed inlet flux:   %.12e mm^3/s", prescribedInletFlowFV1))
print(string.format("Inlet integration rel. error:%.12e", inletIntegrationError))
print(string.format("FV1 nodal-P1 inlet flux:     %.12e mm^3/s", inletFlowFV1))
print(string.format("FV1 outlet flux:             %.12e mm^3/s", outletFlowFV1))
print(string.format("FV1 no-slip Apical flux:     %.12e mm^3/s", apicalFlowFV1))
print(string.format("FV1 flux balance:            %.12e mm^3/s", flowBalanceFV1))
print(string.format("FV1 relative flux balance:   %.12e", relativeFlowBalanceFV1))
print(string.format("Order-5 inlet flux:          %.12e mm^3/s", inletFlowHighOrder))
print(string.format("Order-5 outlet flux:         %.12e mm^3/s", outletFlowHighOrder))
print(string.format("Order-5 no-slip Apical flux: %.12e mm^3/s", apicalFlowHighOrder))
print(string.format("Order-5 flux balance:        %.12e mm^3/s", flowBalanceHighOrder))
print(string.format("Order-5 relative balance:    %.12e",
    relativeFlowBalanceHighOrder))
print(string.format(
    "FV1 PHYSICAL FLUX VALIDATION: %s (relative tolerance %.3e)",
    relativeFlowBalanceFV1 <= fluxBalanceTolerance and "PASS" or "FAIL",
    fluxBalanceTolerance))

print("Stokes refinement control output prefix: " .. outputPrefix)
