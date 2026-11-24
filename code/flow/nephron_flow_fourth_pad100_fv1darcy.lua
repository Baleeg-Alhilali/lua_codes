----------------------------------------------------------
--  FV1 Darcy (comparison; slow solver)
--  True-scale full fourth trace (4th_trace_167.swc, 1 um voxels), 100 um
--  padding mesh in inputs/fourth_trace_full_pad100/.  Derived 2026-11-24 from
--  code/flow/nephron_stokes_darcy_clean.lua: grid/inlet for this mesh,
--  K_mem 7e-12 mm^2 (~9% reabsorption), Outlet-only outflow fix, test
--  options -membranePermeability and -outerPressure, per-node interface
--  flux dump next to the output.  See docs/FOURTH_TRACE_FULL_FLOW_2026-11-24.md.
--  Run (8 ranks):  mpirun -np 8 <UG4>/bin/ugshell -ex code/flow/nephron_flow_fourth_pad100_fv1darcy.lua -numRefs 2
----------------------------------------------------------
--
--  FV1 DARCY TEST (2026-11-24): membrane/Inter Darcy as ConvectionDiffusionFV1,
--  ILU-Uzawa GMG with GMRES(100); Outlet-only outflow fix.
--  TEST COPY (2026-11-24) of code/flow/nephron_stokes_darcy_clean.lua for the
--  true-scale full fourth trace, 100 um padding mesh, K_mem 7e-12 mm^2.  Changes: grid path, inlet centre/direction,
--  output prefix, mean inlet/outlet pressure and DoF statistics printout.
--
--  Steady Stokes--Darcy flow in a projected nephron mesh
--
--  Lumen:       Stokes FV1
--  Membrane:    Darcy FE
--  Interstitium: Darcy FE
--  Interface:   pressure trace, normal traction and BJS slip
--
----------------------------------------------------------

-- Resolve the Model root from this script's location (works on any machine;
-- run from any directory).
local __script = debug.getinfo(1, "S").source
if __script:sub(1, 1) == "@" then __script = __script:sub(2) end
if __script:sub(1, 1) ~= "/" then
    __script = CurrentWorkingDirectory() .. "/" .. __script
end
MODEL_ROOT = assert(__script:match("^(.*)/code/flow/[^/]+$"),
    "Expected this script inside Model/code/flow: " .. __script)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3

numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
gridName = MODEL_ROOT .. "/inputs/fourth_trace_full_pad100/fourth_trace_full_pad100.ugx"
-- Previous mesh: gridName = "runs/third_trace_restored_t3.ugx"
outputName = util.GetParam(
    "-output", MODEL_ROOT .. "/runs/flow_fourth_full_pad100/fv1_r" .. numRefs,
    "VTK output prefix"
)

-- Internal unit system: millimetres and seconds.
-- Velocity: mm/s; kinematic pressure: mm^2/s^2; permeability: mm^2.

-- Inlet geometry and flow rate.  The paper supplies Q, while the 3-D Stokes
-- boundary needs velocity.  For a circular Poiseuille profile,
-- u_mean = Q/(pi R^2) and u_max = 2*u_mean.
-- fourth_trace_full_pad100.ugx (true scale, 1 um voxels, 100 um padding):
-- Inlet on the x-min wall
inletCenterX = 0.360781065
inletCenterY = 0.393057145
inletCenterZ = 0.162855950
inletDirectionX = 1.0
inletDirectionY = 0.0
inletDirectionZ = 0.0
inletRadius = 1.0e-2
inletFlowRate = 5.0e-4                    -- [mm^3/s] = 30 nL/min
meanInletVelocity = inletFlowRate /
    (math.pi * inletRadius * inletRadius)
maximumInletVelocity = 2.0 * meanInletVelocity

-- Previous third_trace_restored_t3.ugx inlet:
-- center = (7.82522082, 11.24569434, 5.583) mm
-- direction = (0.0, 0.0, 1.0), radius = 1.5e-2 mm

-- Fluid and porous-medium parameters
kinematicViscosity = 1.0                  -- [mm^2/s]
membranePermeability = util.GetParamNumber(
    "-membranePermeability", 7.0e-12)     -- [mm^2] (~10% reabsorption, 2026-11-24)
interPermeability = util.GetParamNumber(
    "-interPermeability", 1.0e-11)        -- [mm^2]
bjsPermeability = util.GetParamNumber(
    "-bjsPermeability", 1.0e-13)          -- [mm^2]
bjsAlpha = 0.2                             -- dimensionless
normalImpedance = 1.0e3                    -- [mm/s]
pressureScale = 1.0e4                      -- dimensionless
outerKinematicPressure = util.GetParamNumber(
    "-outerPressure", 6.132812e5)        -- [mm^2/s^2] OuterWall (test option)
outletKinematicPressure = 2.6664474e6      -- [mm^2/s^2]

-- Parallel distribution and validation
firstDistributionLevel = 0
interfaceConservationTolerancePercent = 1.0e-1

--------------------------------------------------------------------------------
-- Initialization
--------------------------------------------------------------------------------

totalStartTime = GetClockS()

InitUG(dim, AlgebraType("CPU", 1))

print("Steady nephron Stokes--Darcy flow")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())
print(string.format(
    "Inlet: Q=%.6e mm^3/s (%.6g nL/min), R=%.6e mm, " ..
    "u_mean=%.6e mm/s, u_max=%.6e mm/s",
    inletFlowRate, inletFlowRate * 6.0e4, inletRadius,
    meanInletVelocity, maximumInletVelocity
))

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
-- Optional extension hook used by the quasi-steady transport driver.  The
-- default flow model is unchanged when this table is absent.
if flowExtraApproximationFunctions ~= nil then
    for _, fct in ipairs(flowExtraApproximationFunctions) do
        approxSpace:add_fct(
            fct.name, "Lagrange", 1, fct.subsets
        )
    end
end
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic("component")

print("Approximation space initialized")

--------------------------------------------------------------------------------
-- Inlet profile
--------------------------------------------------------------------------------

function InletVelocity(x, y, z, t)
    local dx = x - inletCenterX
    local dy = y - inletCenterY
    local dz = z - inletCenterZ
    local axial = dx * inletDirectionX
        + dy * inletDirectionY + dz * inletDirectionZ
    local radialX = dx - axial * inletDirectionX
    local radialY = dy - axial * inletDirectionY
    local radialZ = dz - axial * inletDirectionZ
    local radiusSquared = radialX * radialX
        + radialY * radialY + radialZ * radialZ
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )

    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end

    return maximumInletVelocity * profile * inletDirectionX,
        maximumInletVelocity * profile * inletDirectionY,
        maximumInletVelocity * profile * inletDirectionZ
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
-- FIX TEST (2026-11-24): the outflow condition's continuity term (u.n|face|)
-- on Apical was a mass sink after the pressure-trace merge; Apical keeps the
-- Darcy traction (Neumann), BJS friction and normal impedance only.
stokesOpenBoundary:add("Outlet")

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
-- Keep optional fields fixed during the initial standalone steady-flow solve.
-- The coupled driver builds a new transient domain discretization afterwards.
if flowExtraApproximationFunctions ~= nil then
    steadyExtraFieldBoundary = DirichletBoundary()
    for _, fct in ipairs(flowExtraApproximationFunctions) do
        steadyExtraFieldBoundary:add(0.0, fct.name, fct.subsets)
    end
    domainDisc:add(steadyExtraFieldBoundary)
end
domainDisc:add(OneSideP1Constraints())

darcyPressure = pressureScale * darcyMembraneEq:value()

for _, component in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    interfaceTraction = NeumannBoundary(component[1], "fv1")
    interfaceTraction:add(
        darcyPressure * component[2], "Apical", "Lumen"
    )
    domainDisc:add(interfaceTraction)

    outletTraction = NeumannBoundary(component[1], "fv1")
    outletTraction:add(
        outletKinematicPressure * component[2], "Outlet", "Lumen"
    )
    domainDisc:add(outletTraction)
end

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

-- FV1 Darcy (2026-11-20 two-nephron test): line Gauss-Seidel/BiCGStab did
-- not converge; ILU forward smoother with GMRES(100) did.
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

schurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
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
pressureTrace:set_interface_defect_dump(outputName .. "_interface")
pressureTrace:capture_next_interface_defect()
domainDisc:assemble_defect(interfaceDefect, solution)
interfaceDefectSums = pressureTrace:interface_defect_sums()
stokesInterfaceResidualFlux = interfaceDefectSums[1]
darcyInterfaceResidualFlux = interfaceDefectSums[2]
interfaceFluxMismatch = interfaceDefectSums[3]

--------------------------------------------------------------------------------
-- Results and output
--------------------------------------------------------------------------------

physicalSolution = solution:clone()
physicalScaleFactors = {1.0, 1.0, 1.0, 1.0, pressureScale}
if flowExtraApproximationFunctions ~= nil then
    for _ = 1, #flowExtraApproximationFunctions do
        physicalScaleFactors[#physicalScaleFactors + 1] = 1.0
    end
end
ScaleGF(
    physicalSolution, solution,
    physicalScaleFactors
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

inletPressureMean = Integral(
    physicalSolution, "p_stokes", "Inlet"
) / Integral(1.0, physicalSolution, "Inlet")
outletPressureMean = Integral(
    physicalSolution, "p_stokes", "Outlet"
) / Integral(1.0, physicalSolution, "Outlet")
pressureDrop = inletPressureMean - outletPressureMean

interfaceFluxScale = math.max(
    math.abs(stokesInterfaceResidualFlux),
    math.abs(darcyInterfaceResidualFlux), 1.0e-30
)
interfaceFluxMismatchPercent =
    100.0 * math.abs(interfaceFluxMismatch) / interfaceFluxScale

flowBalance = inletFlow + outletFlow + apicalFlow
flowBalancePercent = 100.0 * math.abs(flowBalance) / math.abs(inletFlow)

print(string.format(
    "Flow [mm^3/s]: inlet=%.6e, outlet=%.6e, apical=%.6e",
    inletFlow, outletFlow, apicalFlow
))
print(string.format(
    "Mean pressure [mm^2/s^2]: inlet=%.12e, outlet=%.12e, drop=%.12e",
    inletPressureMean, outletPressureMean, pressureDrop
))
print(string.format(
    "Lumen balance: %.6e mm^3/s (%.4f%%)",
    flowBalance, flowBalancePercent
))
print(string.format(
    "Assembled interface rows [mm^3/s]: Stokes=%.6e, Darcy=%.6e",
    stokesInterfaceResidualFlux, darcyInterfaceResidualFlux
))
print(string.format(
    "Assembled interface conservation: sum=%.6e mm^3/s (%.6g%%)",
    interfaceFluxMismatch, interfaceFluxMismatchPercent
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

