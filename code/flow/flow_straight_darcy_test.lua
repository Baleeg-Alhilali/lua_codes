--------------------------------------------------------------------------------
-- Darcy-only refinement test
--------------------------------------------------------------------------------
-- Isolates the scalar Darcy FE operator on the same no-bend mesh and with the
-- same permeability used by flow_straight_geometry.lua. A smooth pressure is
-- prescribed on Apical and a constant pressure on OuterWall. This deliberately
-- removes Stokes and PressureTraceConstraint from the system.

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
AssertPluginsLoaded({"ConvectionDiffusion", "neuro_collection", "SuperLU6"})

dim = 3
numRefs = util.GetParamNumber("-numRefs", 0, "number of uniform refinements")
outputPrefix = util.GetParam("-output",
    "/tmp/flow_straight_darcy_ref" .. numRefs, "VTK output prefix")
gridName = "runs/straight_nephron_no_bend.ugx"

viscosity = 1.0                    -- [mm^2/s]
permeabilityM2 = 1.0e-12          -- [m^2]
pressureScale = 1.0e4             -- same stored-pressure scaling as coupled run
density = 1000.0                   -- [kg/m^3]
outerPressurePa = 2.0              -- [Pa]
apicalPressureInletPa = 6.0        -- [Pa]
apicalPressureOutletPa = 1.0       -- [Pa]
axialMin = -0.420                  -- [mm]
axialMax = 0.420                   -- [mm]

InitUG(dim, AlgebraType("CPU", 1))

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {
    "Membrane", "Inter", "Apical", "Basolateral", "OuterWall",
    "MembraneWall", "TripleJunction"
}), "Darcy test mesh is missing a required subset")

balancer.firstDistLvl = numRefs
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, numRefs)
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("p_darcy", "Lagrange", 1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,TripleJunction")
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

local permeability = permeabilityM2 * 1.0e6 -- [mm^2]
local mobility = permeability / viscosity * pressureScale
local outerPressureKinematic =
    outerPressurePa / (density * 1.0e-6) / pressureScale

function smoothApicalPressure(x, y, z, t)
    local s = math.max(0.0, math.min(1.0,
        (z - axialMin) / (axialMax - axialMin)))
    local pressurePa = apicalPressureInletPa
        + s * (apicalPressureOutletPa - apicalPressureInletPa)
    return pressurePa / (density * 1.0e-6) / pressureScale
end

darcyMembrane = ConvectionDiffusionFE("p_darcy", "Membrane")
darcyMembrane:set_mass_scale(0.0)
darcyMembrane:set_diffusion(mobility)
darcyInter = ConvectionDiffusionFE("p_darcy", "Inter")
darcyInter:set_mass_scale(0.0)
darcyInter:set_diffusion(mobility)

pressureBoundary = DirichletBoundary()
pressureBoundary:add("smoothApicalPressure", "p_darcy", "Apical")
pressureBoundary:add(outerPressureKinematic, "p_darcy", "OuterWall")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(darcyMembrane)
domainDisc:add(darcyInter)
domainDisc:add(pressureBoundary)

OrderLex(approxSpace, "+x")
smoother = LineGaussSeidel(approxSpace)
smoother:set_num_steps(1, 1, 1, 1, 1, 1)

matrix = AssembledLinearOperator(domainDisc)
rhs = GridFunction(approxSpace)
solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)
domainDisc:assemble_linear(matrix, rhs)

solver = util.solver.CreateLinearSolver({
    type = "bicgstab",
    precond = {
        type = "gmg", approxSpace = approxSpace,
        discretization = domainDisc, baseLevel = 0,
        baseSolver = "superlu", smoother = smoother,
        cycle = "V", preSmooth = 2, postSmooth = 2, rap = false,
        gatheredBaseSolverIfAmbiguous = true,
        transfer = {type = "std", restrictionDamp = 1.0,
            prolongationDamp = 1.0, enableP1LagrangeOptimization = true}
    },
    convCheck = {type = "standard", iterations = 400,
        absolute = 1.0e-12, reduction = 1.0e-12, verbose = true}
})
solver:init(matrix, solution)
assert(solver:apply(solution, rhs), "Darcy-only GMG solve failed")
solution:set_consistent_storage_type()

defect = rhs:clone()
VecAssign(defect, rhs)
matrix:apply_sub(defect, solution)
print(string.format("Darcy-only true defect: %.12e", VecNorm(defect)))

output = VTKOutput()
output:clear_selection()
output:select_nodal("p_darcy", "darcy_pressure_kinematic")
output:print_subsets(outputPrefix .. "_membrane", solution, "Membrane", 0, 0.0)
output:print_subsets(outputPrefix .. "_inter", solution, "Inter", 0, 0.0)
