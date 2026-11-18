--------------------------------------------------------------------------------
-- Diagnostic summary
--------------------------------------------------------------------------------
-- Isolates refinement, projection, and the scalar Darcy FV1 operator from the
-- Stokes--Darcy interface coupling. A smooth axial pressure is prescribed on
-- Apical and a constant pressure is prescribed on OuterWall. The problem is
-- solved directly, so any coarse/fine pressure oscillation in these outputs
-- must originate in the mesh/projector or Darcy discretization, not GMG,
-- Stokes stabilization, or PressureTraceConstraint.

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
AssertPluginsLoaded({"ConvectionDiffusion", "neuro_collection", "SuperLU6"})

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3
numRefs = util.GetParamNumber("-numRefs", 0, "number of uniform refinements")
gridName = "runs/triplejunction_ogrid40_test.ugx"
outputPrefix = util.GetParam(
    "-output",
    "Results_final/debug_ogrid40_full_darcy_fv1_ref" .. numRefs,
    "VTK output prefix"
)

viscosity = 1.0                    -- [mm^2/s]
membranePermeabilityM2 = 1.0e-9    -- [m^2]
interPermeabilityM2 = 1.0e-8       -- [m^2]
density = 1000.0                   -- [kg/m^3]
outerPressurePa = 4.6e-3           -- [Pa]
apicalPressureInletPa = 5.3e-3     -- [Pa]
apicalPressureOutletPa = 4.8e-3    -- [Pa]
axialMin = 0.2108                  -- coarse-mesh coordinate [mm]
axialMax = 0.3022                  -- coarse-mesh coordinate [mm]

assert(numRefs >= 0, "numRefs must be non-negative")
assert(axialMax > axialMin, "invalid axial interval")

InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {
    "Membrane", "Inter", "Apical", "Basolateral", "OuterWall",
    "MembraneWall", "TripleJunction"
}), "Darcy control mesh is missing a required subset")

balancer.firstDistLvl = 0
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, numRefs)
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(
    "p_darcy", "Lagrange", 1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,TripleJunction"
)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Darcy discretization and boundary conditions
--------------------------------------------------------------------------------

local membranePermeability = membranePermeabilityM2 * 1.0e6 -- [mm^2]
local interPermeability = interPermeabilityM2 * 1.0e6       -- [mm^2]
local membraneMobility = membranePermeability / viscosity
local interMobility = interPermeability / viscosity
local outerPressureKinematic = outerPressurePa / (density * 1.0e-6)

function smoothApicalPressure(x, y, z, t)
    local s = math.max(0.0, math.min(1.0,
        (z - axialMin) / (axialMax - axialMin)))
    local pressurePa = apicalPressureInletPa
        + s * (apicalPressureOutletPa - apicalPressureInletPa)
    return pressurePa / (density * 1.0e-6)
end

darcyMembrane = ConvectionDiffusionFV1("p_darcy", "Membrane")
darcyMembrane:set_mass_scale(0.0)
darcyMembrane:set_diffusion(membraneMobility)

darcyInter = ConvectionDiffusionFV1("p_darcy", "Inter")
darcyInter:set_mass_scale(0.0)
darcyInter:set_diffusion(interMobility)

pressureBoundary = DirichletBoundary()
pressureBoundary:add("smoothApicalPressure", "p_darcy", "Apical")
pressureBoundary:add(outerPressureKinematic, "p_darcy", "OuterWall")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(darcyMembrane)
domainDisc:add(darcyInter)
domainDisc:add(pressureBoundary)

--------------------------------------------------------------------------------
-- Direct solver
--------------------------------------------------------------------------------

matrix = AssembledLinearOperator(domainDisc)
rhs = GridFunction(approxSpace)
solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)
domainDisc:assemble_linear(matrix, rhs)

solver = AgglomeratingSolver(SuperLU())
solver:init(matrix, solution)
assert(solver:apply(solution, rhs), "Darcy control direct solve failed")
solution:set_consistent_storage_type()

defect = rhs:clone()
VecAssign(defect, rhs)
matrix:apply_sub(defect, solution)
print(string.format("Darcy control true defect: %.12e", VecNorm(defect)))

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

output = VTKOutput()
output:clear_selection()
output:select_nodal("p_darcy", "darcy_pressure_kinematic")
output:print_subsets(outputPrefix .. "_membrane", solution, "Membrane", 0, 0.0)
output:print_subsets(outputPrefix .. "_inter", solution, "Inter", 0, 0.0)

print("Darcy refinement control output prefix: " .. outputPrefix)
