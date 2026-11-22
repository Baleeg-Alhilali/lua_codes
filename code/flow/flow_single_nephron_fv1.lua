--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Solves steady three-dimensional flow in one projected nephron. The Lumen
-- uses Stokes FV1, the Membrane and Inter domains use Darcy FV1, and the
-- Apical interface uses scaled pressure trace, normal traction, and BJS slip.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3

numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
gridName = "runs/third_trace_mesh_hex.ugx"
-- Previous mesh: gridName = "runs/third_trace_restored_t3.ugx"
outputName = util.GetParam(
    "-output", "Results_final/single_nephron_fv1_r" .. numRefs,
    "VTK output prefix"
)

-- Internal unit system: millimetres and seconds.
-- Velocity: mm/s; kinematic pressure: mm^2/s^2; permeability: mm^2.

-- Inlet geometry and velocity
inletRadius = 1.0e-2
maximumInletVelocity = 3.0

inletCenter = {
    7.757994361, 11.344866557, 5.588
}
inletDirection = {
    0.0, 0.0, 1.0
}

-- Previous third_trace_restored_t3.ugx inlet:
-- center = (7.82522082, 11.24569434, 5.583) mm
-- direction = (0.0, 0.0, 1.0), radius = 1.5e-2 mm

-- Fluid and porous-medium parameters
kinematicViscosity = 1.0                  -- [mm^2/s]
-- Inherited from the two-nephron FV1 candidate; not recalibrated here.
membranePermeability = 2.35e-12           -- [mm^2], 23.5x baseline
interPermeability = 2.35e-10              -- [mm^2], 23.5x baseline
bjsPermeability = 1.0e-13                 -- [mm^2]
bjsAlpha = 0.2                             -- dimensionless
normalImpedance = 1.0e3                    -- [mm/s]
pressureScale = 1.0e4                      -- dimensionless
outerKinematicPressure = 6.132812e5        -- [mm^2/s^2]
outletKinematicPressure = 2.6664474e6      -- [mm^2/s^2]

-- Parallel distribution and validation
firstDistributionLevel = 0
interfaceConservationTolerancePercent = 1.0e-1
algebraicDefectTolerance = 1.0e-10

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Initializes distributed UG4, refines the projector-backed mesh, and builds
-- the common P1 approximation space for velocity and the two pressures.

totalStartTime = GetClockS()

InitUG(dim, AlgebraType("CPU", 1))

print("Single-nephron FV1 Stokes--Darcy flow")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)

balancer.firstDistLvl = firstDistributionLevel
balancer.parametersParsed = true

print("Refining and distributing...")
balancer.RefineAndRebalanceDomain(dom, numRefs)

print("Domain info:")
print(dom:domain_info():to_string())

allSubsets =
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall," ..
    "MembraneWall,Inlet,Outlet,TripleJunction"

velocitySubsets = {
    "Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"
}

stokesPressureSubsets = {
    "Lumen", "Membrane", "Apical", "Basolateral", "MembraneWall",
    "Inlet", "Outlet", "TripleJunction"
}

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct(
    {"u", "v", "w"}, "Lagrange", 1, velocitySubsets
)
approxSpace:add_fct(
    {"p_stokes"}, "Lagrange", 1, stokesPressureSubsets
)
approxSpace:add_fct(
    "p_darcy_scaled", "Lagrange", 1, allSubsets
)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space initialized")

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Defines the parabolic inlet and registers Stokes FV1, Darcy FV1, pressure
-- trace, traction, BJS, wall, outlet, and inactive-field conditions.

function ParabolicInletVelocity(x, y, z, center, direction)
    local dx = x - center[1]
    local dy = y - center[2]
    local dz = z - center[3]
    local axial = dx * direction[1]
        + dy * direction[2] + dz * direction[3]
    local radialX = dx - axial * direction[1]
    local radialY = dy - axial * direction[2]
    local radialZ = dz - axial * direction[3]
    local radiusSquared = radialX * radialX
        + radialY * radialY + radialZ * radialZ
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )

    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end

    return maximumInletVelocity * profile * direction[1],
        maximumInletVelocity * profile * direction[2],
        maximumInletVelocity * profile * direction[3]
end

function InletVelocity(x, y, z, t)
    return ParabolicInletVelocity(
        x, y, z, inletCenter, inletDirection
    )
end

membraneMobility = membranePermeability / kinematicViscosity
interMobility = interPermeability / kinematicViscosity
scaledMembraneMobility = membraneMobility * pressureScale
scaledInterMobility = interMobility * pressureScale
bjsFriction = bjsAlpha * kinematicViscosity / math.sqrt(bjsPermeability)

outerPressureScaled = outerKinematicPressure / pressureScale

stokesEq = NavierStokesFV1(
    {"u", "v", "w", "p_stokes"}, {"Lumen"}
)
stokesEq:set_stokes(true)
stokesEq:set_laplace(false)
stokesEq:set_kinematic_viscosity(kinematicViscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_upwind("lps")
stokesEq:set_peclet_blend(false)
stokesEq:set_stabilization("fields", "cor")
stokesEq:set_pac_upwind(false)

stokesInlet = NavierStokesInflow(stokesEq)
stokesInlet:add("InletVelocity", "Inlet")

stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
stokesOpenBoundary:add("Outlet,Apical")

stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalImpedance)

darcyMembraneEq = ConvectionDiffusionFV1(
    "p_darcy_scaled", "Membrane"
)
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneMobility)

darcyInterEq = ConvectionDiffusionFV1(
    "p_darcy_scaled", "Inter"
)
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterMobility)

darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(
    outerPressureScaled, "p_darcy_scaled", "OuterWall"
)

inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(
    0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet"
)

inactiveStokes = DirichletBoundary()
inactiveStokes:add(
    0.0, "p_stokes", "Membrane,Basolateral,MembraneWall"
)

pressureTrace = PressureTraceConstraint(
    "p_stokes", "p_darcy_scaled", "Apical,TripleJunction"
)
pressureTrace:set_pressure_scale(pressureScale)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(stokesEq)
domainDisc:add(stokesInlet)
domainDisc:add(stokesOpenBoundary)
domainDisc:add(stokesRim)
domainDisc:add(stokesBJS)
domainDisc:add(darcyMembraneEq)
domainDisc:add(darcyInterEq)
domainDisc:add(pressureTrace)
domainDisc:add(darcyOuterPressure)
domainDisc:add(inactiveDarcy)
domainDisc:add(inactiveStokes)
domainDisc:add(OneSideP1Constraints())

darcyPressure = pressureScale * darcyMembraneEq:value()

for _, component in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    interfaceTraction = NeumannBoundary(component[1], "fv1")
    interfaceTraction:add(
        darcyPressure * component[2],
        "Apical", "Lumen"
    )
    domainDisc:add(interfaceTraction)

    outletTraction = NeumannBoundary(component[1], "fv1")
    outletTraction:add(
        outletKinematicPressure * component[2],
        "Outlet", "Lumen"
    )
    domainDisc:add(outletTraction)
end

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

mainInverse = ILU()
mainInverse:set_damp(1.0)
mainInverse:set_inversion_eps(1.0e-30)

pressureInverse = ILU()
pressureInverse:set_damp(1.0)
pressureInverse:set_inversion_eps(1.0e-30)

uzawa = UzawaBase({"p_stokes"})
uzawa:set_forward_iter(mainInverse)
uzawa:set_schur_iter(pressureInverse)
uzawa:set_backward_iter(mainInverse)

schurMassEq = ConvectionDiffusionFV1(
    "p_stokes", "Lumen"
)
schurMassEq:set_mass_scale(1.0)
schurMassEq:set_diffusion(0.0)

schurMassDisc = DomainDiscretization(approxSpace)
schurMassDisc:add(schurMassEq)
schurMassDisc:add(OneSideP1Constraints())

schurMassOp = AssembledLinearOperator(schurMassDisc)
uzawa:set_schur_operator_update(
    schurMassOp, 1.0 / kinematicViscosity
)

solverDesc =
{
    type = "gmres", -- GMRES(100) is constructed directly below
    precond =
    {
        type = "gmg", -- verified scalable multilevel preconditioner
        approxSpace = approxSpace,
        discretization = domainDisc,
        rap = true, -- verified Galerkin coarse operators for this coupling
        smoother = uzawa, -- verified Stokes pressure-block smoother
        cycle = "V",
        preSmooth = 2,
        postSmooth = 2,
        baseSolver = "superlu", -- coarse grid only; not the refined solver
        baseLevel = 0,
        gatheredBaseSolverIfAmbiguous = true,
        transfer =
        {
            type = "std",
            restrictionDamp = 1.0,
            prolongationDamp = 1.0,
            enableP1LagrangeOptimization = true
        }
    },
    convCheck =
    {
        type = "standard",
        iterations = 400,
        absolute = 1.0e-14,
        reduction = 1.0e-8,
        verbose = true
    }
}

util.solver.defaults.approxSpace = approxSpace
solver = GMRES(100)
solver:set_preconditioner(
    util.solver.CreatePreconditioner(solverDesc.precond)
)
solver:set_convergence_check(
    util.solver.CreateConvCheck(solverDesc.convCheck)
)

matrix = AssembledLinearOperator(domainDisc)

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- This is a steady linear problem: no time integrator is used, and the
-- algebraic solution starts from zero before boundary adjustment.

solution = GridFunction(approxSpace)
rightHandSide = solution:clone()
solution:set(0.0)

domainDisc:adjust_solution(solution)
solution:set_consistent_storage_type()
domainDisc:assemble_linear(matrix, rightHandSide)

print("Solving...")
solveStartTime = GetClockS()

solver:init(matrix, solution)
if not solver:apply(solution, rightHandSide) then
    print("Solver failed")
    exit()
end

solveTime = GetClockS() - solveStartTime
solution:set_consistent_storage_type()

-- Independently recompute the true assembled residual and its five components.
trueDefect = rightHandSide:clone()
VecAssign(trueDefect, rightHandSide)
matrix:apply_sub(trueDefect, solution)
trueDefectNorm = VecNorm(trueDefect)
componentNames = {"u", "v", "w", "p_stokes", "p_darcy_scaled"}
componentDefectNorms = {}
for component = 1, #componentNames do
    local factors = {0.0, 0.0, 0.0, 0.0, 0.0}
    factors[component] = 1.0
    local selectedDefect = trueDefect:clone()
    ScaleGF(selectedDefect, trueDefect, factors)
    componentDefectNorms[component] = VecNorm(selectedDefect)
end
print(string.format("Recomputed true defect: %.12e", trueDefectNorm))
for component = 1, #componentNames do
    print(string.format(
        "  %s defect: %.12e",
        componentNames[component], componentDefectNorms[component]
    ))
end
if trueDefectNorm > algebraicDefectTolerance then
    print(string.format(
        "Algebraic validation failed: %.6e > %.6e",
        trueDefectNorm, algebraicDefectTolerance
    ))
    exit()
end

-- Capture the separate interface rows before PressureTraceConstraint merges
-- them. This measures the exact algebraic normal-flux conservation equation.
interfaceDefect = solution:clone()
pressureTrace:capture_next_interface_defect()
domainDisc:assemble_defect(interfaceDefect, solution)
interfaceDefectSums = pressureTrace:interface_defect_sums()

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- Converts stored Darcy pressure to physical kinematic pressure, reports
-- algebraic and conservation diagnostics, and writes all three flow domains.

physicalSolution = solution:clone()
ScaleGF(
    physicalSolution, solution,
    {1.0, 1.0, 1.0, 1.0, pressureScale}
)
physicalSolution:set_consistent_storage_type()

stokesVelocity = stokesEq:velocity()

darcyMembraneVelocity = DarcyVelocityLinker()
darcyMembraneVelocity:set_permeability(membranePermeability)
darcyMembraneVelocity:set_viscosity(kinematicViscosity)
darcyMembraneVelocity:set_density(1.0)
darcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
darcyMembraneVelocity:set_pressure_gradient(darcyMembraneEq:gradient())

darcyInterVelocity = DarcyVelocityLinker()
darcyInterVelocity:set_permeability(interPermeability)
darcyInterVelocity:set_viscosity(kinematicViscosity)
darcyInterVelocity:set_density(1.0)
darcyInterVelocity:set_gravity(ConstUserVector(0.0))
darcyInterVelocity:set_pressure_gradient(darcyInterEq:gradient())

analyticInletFlow = 0.5 * math.pi * inletRadius * inletRadius
    * maximumInletVelocity

inletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalSolution, "Inlet", "Lumen"
)
outletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalSolution, "Outlet", "Lumen"
)
apicalFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalSolution, "Apical", "Lumen"
)
darcyApicalFlow = IntegralNormalComponentOnManifold(
    darcyMembraneVelocity, physicalSolution, "Apical", "Membrane"
)

flowBalance = inletFlow + outletFlow + apicalFlow
flowBalancePercent = 100.0 * math.abs(flowBalance) / math.abs(inletFlow)
outletPercent = 100.0 * outletFlow / analyticInletFlow
membranePercent = 100.0 - outletPercent

interfaceFluxScale = math.max(
    math.abs(interfaceDefectSums[1]),
    math.abs(interfaceDefectSums[2]), 1.0e-30
)
interfaceFluxMismatchPercent =
    100.0 * math.abs(interfaceDefectSums[3]) / interfaceFluxScale

print(string.format(
    "Flow [mm^3/s]: inlet=%.6e, outlet=%.6e, apical=%.6e",
    inletFlow, outletFlow, apicalFlow
))
print(string.format(
    "Lumen balance: %.6e mm^3/s (%.4f%%)",
    flowBalance, flowBalancePercent
))
print(string.format(
    "Analytic split: outlet=%.3f%%, membrane=%.3f%%",
    outletPercent, membranePercent
))
print(string.format(
    "Assembled interface rows [mm^3/s]: Stokes=%.6e, Darcy=%.6e",
    interfaceDefectSums[1], interfaceDefectSums[2]
))
print(string.format(
    "Assembled interface conservation: sum=%.6e mm^3/s (%.6g%%)",
    interfaceDefectSums[3], interfaceFluxMismatchPercent
))
print(string.format(
    "Reconstructed Apical fluxes [mm^3/s]: Stokes=%.6e, Darcy=%.6e",
    apicalFlow, darcyApicalFlow
))

if interfaceFluxMismatchPercent > interfaceConservationTolerancePercent then
    print(string.format(
        "Interface conservation failed: %.6g%% > %.6g%%",
        interfaceFluxMismatchPercent,
        interfaceConservationTolerancePercent
    ))
    exit()
end

lumenOutput = VTKOutput()
lumenOutput:clear_selection()
lumenOutput:select_nodal({"u", "v", "w"}, "velocity_mm_s")
lumenOutput:select_nodal("p_stokes", "pressure_kinematic")
lumenOutput:print_subsets(
    outputName .. "_lumen", physicalSolution, "Lumen", 0, 0.0
)

function VectorComponents(vectorData)
    local components = {}
    for component = 0, 2 do
        components[component + 1] = UserVectorEntryAdapter3d()
        components[component + 1]:set_vector(vectorData, component)
    end
    return components
end

function WriteDarcyOutput(suffix, subset, velocity)
    local components = VectorComponents(velocity)
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal(
        "p_darcy_scaled", "pressure_kinematic"
    )
    output:select_element(velocity, "velocity_mm_s")
    output:select_element(components[1], "u")
    output:select_element(components[2], "v")
    output:select_element(components[3], "w")
    output:print_subsets(
        outputName .. suffix, physicalSolution, subset, 0, 0.0
    )
end

WriteDarcyOutput(
    "_membrane", "Membrane", darcyMembraneVelocity
)
WriteDarcyOutput(
    "_inter", "Inter", darcyInterVelocity
)

totalTime = GetClockS() - totalStartTime

print(string.format("Solver time: %.3f s", solveTime))
print(string.format("Total time:  %.3f s", totalTime))
print("Done")




