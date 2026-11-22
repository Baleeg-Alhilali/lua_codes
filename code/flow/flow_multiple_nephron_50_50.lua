----------------------------------------------------------
--
--  Calibrated 50/50 Stokes--Darcy flow in two projected nephrons
--
--  Lumen:       Stokes FV1
--  Membrane:    Darcy FE
--  Interstitium: Darcy FE
--  Interface:   pressure trace, normal traction and BJS slip
--
----------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3

numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
gridName = "runs/third_fourth_nephrons_membranewall_fixed.ugx"
outputName = util.GetParam(
    "-output", "Results_final/multi_nephron_50_50_r" .. numRefs,
    "VTK output prefix"
)

-- Internal unit system: millimetres and seconds.
-- Velocity: mm/s; kinematic pressure: mm^2/s^2; permeability: mm^2.

-- Inlet geometry and velocity
inletRadius = 1.5e-2
maximumInletVelocity = 3.0

inlet1Center = {
    7.473194986785253, 12.454510944000001, 5.6048934974042073
}
inlet1Direction = {
    -1.038317342773153e-5, -0.9999999999425444, 2.6647391101191372e-6
}

inlet2Center = {
    12.207611731404208, 9.941724906285252, 6.7059727200000001
}
inlet2Direction = {
    3.7090925200442324e-11, -7.4346346465763484e-12, -1.0
}

-- Fluid and porous-medium parameters
kinematicViscosity = 1.0                  -- [mm^2/s]
membranePermeability = 2.35e-12           -- [mm^2], calibrated 23.5x baseline
interPermeability = 2.35e-10              -- [mm^2], calibrated 23.5x baseline
bjsPermeability = 1.0e-13                 -- [mm^2]
bjsAlpha = 0.2                             -- dimensionless
normalImpedance = 1.0e3                    -- [mm/s]
pressureScale = 1.0e4                      -- dimensionless
outerKinematicPressure = 6.132812e5        -- [mm^2/s^2]
outletKinematicPressure = 2.6664474e6      -- [mm^2/s^2]

-- Parallel distribution and validation
firstDistributionLevel = 0
interfaceConservationTolerancePercent = 1.0e-1

--------------------------------------------------------------------------------
-- Initialization
--------------------------------------------------------------------------------

totalStartTime = GetClockS()

InitUG(dim, AlgebraType("CPU", 1))

print("Calibrated 50/50 two-nephron Stokes--Darcy flow")
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
    "Lumen_1,Lumen_2,Membrane_1,Membrane_2,Inter," ..
    "Apical_1,Apical_2,Basolateral_1,Basolateral_2,OuterWall," ..
    "MembraneWall_1,MembraneWall_2,Inlet_1,Inlet_2," ..
    "Outlet_1,Outlet_2,TripleJunction_1,TripleJunction_2"

velocitySubsets = {
    "Lumen_1", "Lumen_2", "Apical_1", "Apical_2",
    "Inlet_1", "Inlet_2", "Outlet_1", "Outlet_2",
    "TripleJunction_1", "TripleJunction_2"
}

stokesPressureSubsets = {
    "Lumen_1", "Lumen_2", "Membrane_1", "Membrane_2",
    "Apical_1", "Apical_2", "Basolateral_1", "Basolateral_2",
    "MembraneWall_1", "MembraneWall_2", "Inlet_1", "Inlet_2",
    "Outlet_1", "Outlet_2", "TripleJunction_1", "TripleJunction_2"
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
-- Inlet profile
--------------------------------------------------------------------------------

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

function InletVelocity1(x, y, z, t)
    return ParabolicInletVelocity(
        x, y, z, inlet1Center, inlet1Direction
    )
end

function InletVelocity2(x, y, z, t)
    return ParabolicInletVelocity(
        x, y, z, inlet2Center, inlet2Direction
    )
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------

membraneMobility = membranePermeability / kinematicViscosity
interMobility = interPermeability / kinematicViscosity
scaledMembraneMobility = membraneMobility * pressureScale
scaledInterMobility = interMobility * pressureScale
bjsFriction = bjsAlpha * kinematicViscosity / math.sqrt(bjsPermeability)

outerPressureScaled = outerKinematicPressure / pressureScale

stokesEq = NavierStokesFV1(
    {"u", "v", "w", "p_stokes"}, {"Lumen_1", "Lumen_2"}
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
stokesInlet:add("InletVelocity1", "Inlet_1")
stokesInlet:add("InletVelocity2", "Inlet_2")

stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
stokesOpenBoundary:add("Outlet_1,Outlet_2,Apical_1,Apical_2")

stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction_1,TripleJunction_2")

stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical_1,Apical_2")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalImpedance)

darcyMembraneEq = ConvectionDiffusionFE(
    "p_darcy_scaled", "Membrane_1,Membrane_2"
)
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneMobility)

darcyInterEq = ConvectionDiffusionFE(
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
    0.0, "p_darcy_scaled",
    "Lumen_1,Lumen_2,Inlet_1,Inlet_2,Outlet_1,Outlet_2"
)

inactiveStokes = DirichletBoundary()
inactiveStokes:add(
    0.0, "p_stokes",
    "Membrane_1,Membrane_2,Basolateral_1,Basolateral_2," ..
    "MembraneWall_1,MembraneWall_2"
)

pressureTraceConstraints = {}
for nephron = 1, 2 do
    local suffix = "_" .. nephron
    pressureTraceConstraints[nephron] = PressureTraceConstraint(
        "p_stokes", "p_darcy_scaled",
        "Apical" .. suffix .. ",TripleJunction" .. suffix
    )
    pressureTraceConstraints[nephron]:set_pressure_scale(pressureScale)
end

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(stokesEq)
domainDisc:add(stokesInlet)
domainDisc:add(stokesOpenBoundary)
domainDisc:add(stokesRim)
domainDisc:add(stokesBJS)
domainDisc:add(darcyMembraneEq)
domainDisc:add(darcyInterEq)
for nephron = 1, 2 do
    domainDisc:add(pressureTraceConstraints[nephron])
end
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
        "Apical_1,Apical_2", "Lumen_1,Lumen_2"
    )
    domainDisc:add(interfaceTraction)

    outletTraction = NeumannBoundary(component[1], "fv1")
    outletTraction:add(
        outletKinematicPressure * component[2],
        "Outlet_1,Outlet_2", "Lumen_1,Lumen_2"
    )
    domainDisc:add(outletTraction)
end

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

mainInverse = LineGaussSeidel(approxSpace)
mainInverse:set_num_steps(1, 1, 1, 1, 1, 1)

pressureInverse = ILU()
pressureInverse:set_damp(1.0)
pressureInverse:set_inversion_eps(1.0e-30)

uzawa = UzawaBase({"p_stokes"})
uzawa:set_forward_iter(mainInverse)
uzawa:set_schur_iter(pressureInverse)
uzawa:set_backward_iter(mainInverse)

schurMassEq = ConvectionDiffusionFV1(
    "p_stokes", "Lumen_1,Lumen_2"
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
    type = "bicgstab", -- verified linear outer solver
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
solver = util.solver.CreateLinearSolver(solverDesc)

matrix = AssembledLinearOperator(domainDisc)
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

-- Capture the separate interface rows before PressureTraceConstraint merges
-- them. This measures the exact algebraic normal-flux conservation equation.
interfaceDefect = solution:clone()
for nephron = 1, 2 do
    pressureTraceConstraints[nephron]:capture_next_interface_defect()
end
domainDisc:assemble_defect(interfaceDefect, solution)
interfaceDefectSums = {}
for nephron = 1, 2 do
    interfaceDefectSums[nephron] =
        pressureTraceConstraints[nephron]:interface_defect_sums()
end

--------------------------------------------------------------------------------
-- Results and output
--------------------------------------------------------------------------------

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

nephronFlows = {}
for nephron = 1, 2 do
    local suffix = "_" .. nephron
    local lumen = "Lumen" .. suffix
    local membrane = "Membrane" .. suffix

    nephronFlows[nephron] = {}
    nephronFlows[nephron].inlet = IntegralNormalComponentOnManifold(
        stokesVelocity, physicalSolution, "Inlet" .. suffix, lumen
    )
    nephronFlows[nephron].outlet = IntegralNormalComponentOnManifold(
        stokesVelocity, physicalSolution, "Outlet" .. suffix, lumen
    )
    nephronFlows[nephron].apical = IntegralNormalComponentOnManifold(
        stokesVelocity, physicalSolution, "Apical" .. suffix, lumen
    )
    nephronFlows[nephron].darcyApical =
        IntegralNormalComponentOnManifold(
            darcyMembraneVelocity, physicalSolution,
            "Apical" .. suffix, membrane
        )
    nephronFlows[nephron].balance =
        nephronFlows[nephron].inlet
        + nephronFlows[nephron].outlet
        + nephronFlows[nephron].apical
    nephronFlows[nephron].balancePercent =
        100.0 * math.abs(nephronFlows[nephron].balance)
        / math.abs(nephronFlows[nephron].inlet)
    nephronFlows[nephron].outletPercent =
        100.0 * nephronFlows[nephron].outlet / analyticInletFlow
    nephronFlows[nephron].membranePercent =
        100.0 - nephronFlows[nephron].outletPercent
end

maximumInterfaceFluxMismatchPercent = 0.0
for nephron = 1, 2 do
    local sums = interfaceDefectSums[nephron]
    local fluxScale = math.max(
        math.abs(sums[1]), math.abs(sums[2]), 1.0e-30
    )
    sums.mismatchPercent = 100.0 * math.abs(sums[3]) / fluxScale
    maximumInterfaceFluxMismatchPercent = math.max(
        maximumInterfaceFluxMismatchPercent, sums.mismatchPercent
    )
end

for nephron = 1, 2 do
    local flow = nephronFlows[nephron]
    print(string.format(
        "Nephron %d flow [mm^3/s]: inlet=%.6e, outlet=%.6e, apical=%.6e",
        nephron, flow.inlet, flow.outlet, flow.apical
    ))
    print(string.format(
        "Nephron %d lumen balance: %.6e mm^3/s (%.4f%%)",
        nephron, flow.balance, flow.balancePercent
    ))
    print(string.format(
        "Nephron %d analytic split: outlet=%.3f%%, membrane=%.3f%%",
        nephron, flow.outletPercent, flow.membranePercent
    ))
end
for nephron = 1, 2 do
    local sums = interfaceDefectSums[nephron]
    print(string.format(
        "Nephron %d assembled interface rows [mm^3/s]: " ..
        "Stokes=%.6e, Darcy=%.6e",
        nephron, sums[1], sums[2]
    ))
    print(string.format(
        "Nephron %d assembled interface conservation: " ..
        "sum=%.6e mm^3/s (%.6g%%)",
        nephron, sums[3], sums.mismatchPercent
    ))
end
for nephron = 1, 2 do
    local flow = nephronFlows[nephron]
    print(string.format(
        "Nephron %d reconstructed Apical fluxes [mm^3/s]: " ..
        "Stokes=%.6e, Darcy=%.6e",
        nephron, flow.apical, flow.darcyApical
    ))
end

if maximumInterfaceFluxMismatchPercent
        > interfaceConservationTolerancePercent then
    print(string.format(
        "Interface conservation failed: %.6g%% > %.6g%%",
        maximumInterfaceFluxMismatchPercent,
        interfaceConservationTolerancePercent
    ))
    exit()
end

lumenOutput = VTKOutput()
lumenOutput:clear_selection()
lumenOutput:select_nodal({"u", "v", "w"}, "velocity_mm_s")
lumenOutput:select_nodal("p_stokes", "pressure_kinematic")
lumenOutput:print_subsets(
    outputName .. "_lumen", physicalSolution,
    "Lumen_1,Lumen_2", 0, 0.0
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
    "_membrane", "Membrane_1,Membrane_2", darcyMembraneVelocity
)
WriteDarcyOutput(
    "_inter", "Inter", darcyInterVelocity
)

totalTime = GetClockS() - totalStartTime

print(string.format("Solver time: %.3f s", solveTime))
print(string.format("Total time:  %.3f s", totalTime))
print("Done")
