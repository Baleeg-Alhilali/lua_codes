--------------------------------------------------------------------------------
-- Diagnostic summary
--------------------------------------------------------------------------------
-- Tests Darcy FV1 on the conforming ProMesh Long_tube_V3 grid.  Membrane and
-- Inter use the same permeability, a smooth axial pressure is prescribed on
-- Apical, and a constant pressure is prescribed on OuterWall.  The direct
-- solve excludes Stokes coupling and pressure scaling so that refined-grid
-- pressure patches and Darcy multigrid behavior can be tested independently.

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
AssertPluginsLoaded({"ConvectionDiffusion", "SuperLU6"})

local totalStart = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3
numRefs = util.GetParamNumber("-numRefs", 0, "number of uniform refinements")
gridName = "runs/Long_tube_V3.ugx"
outputPrefix = util.GetParam(
    "-output",
    "Results_final/Long_tube_V3_darcy_fv1_control_ref" .. numRefs,
    "VTK output prefix"
)

-- Long_tube_V3 coordinates are in micrometres: radius 10 um, length 1000 um.
-- The scalar pressure test is invariant to a consistent length-unit choice;
-- only the shape of the pressure field is being diagnosed here.
viscosity = 1.0
permeability = 1.0e-6
outerPressure = 4600.0
apicalPressureInlet = 5300.0
apicalPressureOutlet = 4800.0
axialMin = 0.0
axialMax = 1000.0

assert(numRefs >= 0, "numRefs must be non-negative")
InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------

dom = Domain()
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {
    "Membrane", "Inter", "Apical", "Basolateral", "OuterWall",
    "MembraneWall", "TripleJunction"
}), "Long_tube_V3 is missing a required subset")

-- Keep the coarse hierarchy undistributed.  Distributing level 0 caused the
-- prism hierarchy's parallel coarse-grid correction to be ineffective.
balancer.firstDistLvl = numRefs
balancer.parametersParsed = true
local refineStart = GetClockS()
balancer.RefineAndRebalanceDomain(dom, numRefs)
print(string.format("Mesh load/refinement time: %.6f s", GetClockS() - refineStart))
print("Mesh statistics:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(
    "p_darcy", "Lagrange", 1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,TripleJunction"
)
approxSpace:init_levels()
approxSpace:init_top_surface()
print("Approximation-space statistics:")
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Darcy FV1 discretization and boundary conditions
--------------------------------------------------------------------------------

function smoothApicalPressure(x, y, z, t)
    local s = math.max(0.0, math.min(1.0,
        (z - axialMin) / (axialMax - axialMin)))
    return apicalPressureInlet
        + s * (apicalPressureOutlet - apicalPressureInlet)
end

darcyMembrane = ConvectionDiffusionFV1("p_darcy", "Membrane")
darcyMembrane:set_mass_scale(0.0)
darcyMembrane:set_diffusion(permeability / viscosity)

darcyInter = ConvectionDiffusionFV1("p_darcy", "Inter")
darcyInter:set_mass_scale(0.0)
darcyInter:set_diffusion(permeability / viscosity)

pressureBoundary = DirichletBoundary()
pressureBoundary:add("smoothApicalPressure", "p_darcy", "Apical")
pressureBoundary:add(outerPressure, "p_darcy", "OuterWall")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(darcyMembrane)
domainDisc:add(darcyInter)
domainDisc:add(pressureBoundary)

--------------------------------------------------------------------------------
-- Darcy-only GMG diagnostic solve
--------------------------------------------------------------------------------

matrix = AssembledLinearOperator(domainDisc)
rhs = GridFunction(approxSpace)
solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)

local assemblyStart = GetClockS()
domainDisc:assemble_linear(matrix, rhs)
print(string.format("Assembly time: %.6f s", GetClockS() - assemblyStart))

local solveStart = GetClockS()
local smoother = ElementGaussSeidel(1.0, "element")

local solverDescription = {
    type = "bicgstab",
    precond = {
        type = "gmg",
        approxSpace = approxSpace,
        discretization = domainDisc,
        baseLevel = 0,
        baseSolver = "superlu",
        smoother = smoother,
        cycle = "V",
        preSmooth = 2,
        postSmooth = 2,
        rap = false,
        gatheredBaseSolverIfAmbiguous = true,
        transfer = {
            type = "std",
            restrictionDamp = 1.0,
            prolongationDamp = 1.0,
            enableP1LagrangeOptimization = true
        }
    },
    convCheck = {
        type = "standard",
        iterations = 400,
        absolute = 1.0e-12,
        reduction = 1.0e-10,
        verbose = true
    }
}
util.solver.defaults.approxSpace = approxSpace
solver = util.solver.CreateLinearSolver(solverDescription)
solver:init(matrix, solution)
local solverReportedConvergence = solver:apply(solution, rhs)
solution:set_consistent_storage_type()
print("Solver reported convergence: " .. tostring(solverReportedConvergence))
print(string.format("Solve time: %.6f s", GetClockS() - solveStart))

defect = rhs:clone()
VecAssign(defect, rhs)
matrix:apply_sub(defect, solution)
print(string.format("True absolute defect: %.12e", VecNorm(defect)))

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

output = VTKOutput()
output:clear_selection()
output:select_nodal("p_darcy", "darcy_pressure")
output:print_subsets(outputPrefix .. "_membrane", solution, "Membrane", 0, 0.0)
output:print_subsets(outputPrefix .. "_inter", solution, "Inter", 0, 0.0)

print("Darcy FV1 control output prefix: " .. outputPrefix)
print(string.format("Total runtime: %.6f s", GetClockS() - totalStart))
