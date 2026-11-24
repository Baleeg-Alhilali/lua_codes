----------------------------------------------------------
--
--  Standalone quasi-steady nephron flow and multispecies transport
--
--  Lumen:       Stokes FV1
--  Membrane:    Darcy FE
--  Interstitium: Darcy FE
--  Interface:   pressure trace, normal traction and BJS slip
--
--  Transport: transient Na+, K+, Cl- and urea advection--diffusion
--  Apical/Basolateral: published passive leak pathways only
--  No Na/K pump and no solvent-drag flux in this model stage.
--
----------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3

numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")
gridName = "/Users/alhilaba/Documents/Promesh_swc/Model/runs/third_trace_third_hex_sim_ready.ugx"
-- Previous mesh: gridName = "runs/third_trace_restored_t3.ugx"
outputName = util.GetParam(
    "-output", "Results_final/nephron_flow_transport_r" .. numRefs,
    "VTK output prefix"
)

local function EnsureOutputDirectory(prefix)
    local directory = string.match(prefix, "^(.*)/[^/]*$")
    if directory == nil or directory == "" then return end

    local current = string.sub(directory, 1, 1) == "/" and "/" or ""
    for part in string.gmatch(directory, "[^/]+") do
        if current == "/" then
            current = current .. part
        elseif current == "" then
            current = part
        else
            current = current .. "/" .. part
        end
        if not DirectoryExists(current) then
            assert(CreateDirectory(current) or DirectoryExists(current),
                   "Could not create output directory: " .. current)
        end
    end
end

-- Internal unit system: millimetres and seconds.
-- Velocity: mm/s; kinematic pressure: mm^2/s^2; permeability: mm^2.

-- Inlet geometry and flow rate.  The paper supplies Q, while the 3-D Stokes
-- boundary needs velocity.  For a circular Poiseuille profile,
-- u_mean = Q/(pi R^2) and u_max = 2*u_mean.
inletCenterX = 7.74237
inletCenterY = 11.3605
inletCenterZ = 5.538
inletDirectionX = 0.0
inletDirectionY = 0.0
inletDirectionZ = 1.0
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
    "-membranePermeability", 1.0e-13)     -- [mm^2]
interPermeability = util.GetParamNumber(
    "-interPermeability", 1.0e-11)        -- [mm^2]
bjsPermeability = util.GetParamNumber(
    "-bjsPermeability", 1.0e-13)          -- [mm^2]
bjsAlpha = 0.2                             -- dimensionless
normalImpedance = 1.0e3                    -- [mm/s]
pressureScale = 1.0e4                      -- dimensionless
outerKinematicPressure = 6.132812e5        -- [mm^2/s^2]
outletKinematicPressure = 2.6664474e6      -- [mm^2/s^2]

-- Parallel distribution and validation
firstDistributionLevel = 0
interfaceConservationTolerancePercent = 1.0e-1

-- Transient transport run control
endTime = util.GetParamNumber("-endTime", 100.0)
dt = util.GetParamNumber("-dt", 1.0)
outputEverySteps = util.GetParamNumber("-outputEverySteps", 1)

-- Concentrations [mM] and free-solution diffusion [mm^2/s]. The validated
-- urea formulation uses the same aqueous diffusion coefficient in all three
-- resolved volumes; this leak-only extension follows that formulation rather
-- than introducing an unsupported 1000-fold intracellular hindrance factor.
--
-- Weinstein's h_i,alpha*A coefficients are already normalized by smooth
-- epithelial reference area. The 3-D MembraneTransportFV1 operator integrates
-- flux density over the actual smooth face area, so no radius or 2*pi*r factor
-- is applied here. Apical uses MI. Basolateral combines the parallel IE and IS
-- cell-to-interspace pathways. Zero paper coefficients remain zero; whole-
-- epithelium effective permeabilities are not duplicated on both faces.
soluteParameters = {
    na = {
        label = "Na", diffusion = 1.33e-3, inlet = 140.0,
        cellInitial = 19.6, interInitial = 140.3,
        apicalPermeability = 0.0,
        basolateralPermeability = 1.40e-6,
    },
    k = {
        label = "K", diffusion = 1.96e-3, inlet = 4.9,
        cellInitial = 138.1, interInitial = 4.66,
        apicalPermeability = 9.00e-5,
        basolateralPermeability = 7.40e-4,
    },
    cl = {
        label = "Cl", diffusion = 2.03e-3, inlet = 113.2,
        cellInitial = 16.3, interInitial = 112.0,
        apicalPermeability = 0.0,
        basolateralPermeability = 0.0,
    },
    urea = {
        label = "Urea", diffusion = 1.67e-3, inlet = 5.0,
        cellInitial = 4.96, interInitial = 4.91,
        apicalPermeability = 3.78e-4,
        basolateralPermeability = 3.70e-4,
    },
}
soluteOrder = {"na", "k", "cl", "urea"}
compartmentOrder = {"lumen", "membrane", "inter"}

membraneTransportDensity = 1.0

assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0 and dt <= endTime, "-dt must be in (0, endTime]")
assert(outputEverySteps >= 1 and
       outputEverySteps == math.floor(outputEverySteps),
       "-outputEverySteps must be a positive integer")

--------------------------------------------------------------------------------
-- Initialization
--------------------------------------------------------------------------------

totalStartTime = GetClockS()

InitUG(dim, AlgebraType("CPU", 1))

if ProcRank() == 0 then EnsureOutputDirectory(outputName) end
SynchronizeProcesses()

print("Standalone quasi-steady nephron flow and multispecies transport")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())
print(string.format(
    "Inlet: Q=%.6e mm^3/s (%.6g nL/min), R=%.6e mm, " ..
    "u_mean=%.6e mm/s, u_max=%.6e mm/s",
    inletFlowRate, inletFlowRate * 6.0e4, inletRadius,
    meanInletVelocity, maximumInletVelocity
))
print(string.format("Transport interval: 0 to %.6g s, dt=%.6g s",
                    endTime, dt))
for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    print(string.format(
        "%s: inlet=%.6g cell0=%.6g inter0=%.6g mM; " ..
        "D=%.6g mm^2/s; " ..
        "Papical=%.6g Pbaso=%.6g mm/s",
        p.label, p.inlet, p.cellInitial, p.interInitial,
        p.diffusion,
        p.apicalPermeability, p.basolateralPermeability))
end

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
transportFieldSubsets = {
    lumen = "Lumen,Apical,Inlet,Outlet,TripleJunction",
    membrane = "Membrane,Apical,Basolateral,MembraneWall,TripleJunction",
    inter = "Inter,Basolateral,OuterWall",
}
transportFields = {}
for _, soluteKey in ipairs(soluteOrder) do
    transportFields[soluteKey] = {}
    for _, compartment in ipairs(compartmentOrder) do
        local field = soluteKey .. "_" .. compartment
        transportFields[soluteKey][compartment] = field
        approxSpace:add_fct(
            field, "Lagrange", 1, transportFieldSubsets[compartment])
    end
end
approxSpace:init_levels()
approxSpace:init_top_surface()

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
-- Outlet only (2026-11-24): on Apical, the continuity term (u.n |face|) of
-- this condition was a mass sink after the PressureTraceConstraint merge,
-- because the lumen and membrane control volumes around an interface node
-- form one union CV whose Apical face is internal; the removed fluid never
-- entered the membrane (4.1% of Q at r2 with K_mem 7e-12 mm^2).  Its
-- tangential viscous flux also double-counted with the BJS friction.  Apical
-- keeps the Darcy traction (Neumann), BJS friction and normal impedance.
stokesOpenBoundary:add("Outlet")

stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalImpedance)

darcyMembraneEq = ConvectionDiffusionFE(
    "p_darcy_scaled", "Membrane"
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
-- Fix every concentration field only during the initial steady-flow solve.
-- The transient coupled discretization below replaces this boundary object.
steadyTransportBoundary = DirichletBoundary()
for _, soluteKey in ipairs(soluteOrder) do
    for _, compartment in ipairs(compartmentOrder) do
        steadyTransportBoundary:add(
            0.0, transportFields[soluteKey][compartment],
            transportFieldSubsets[compartment])
    end
end
domainDisc:add(steadyTransportBoundary)
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

mainInverse = LineGaussSeidel(approxSpace)
mainInverse:set_num_steps(1, 1, 1, 1, 1, 1)

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
for _ = 1, #soluteOrder * #compartmentOrder do
    physicalScaleFactors[#physicalScaleFactors + 1] = 1.0
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
print("Steady flow stage complete")

--------------------------------------------------------------------------------
-- Quasi-steady transient transport discretization
--------------------------------------------------------------------------------

-- This follows the validated nephron_flow_urea_transient.lua formulation.
-- Flow and all concentrations share one approximation space, so UG4 can
-- evaluate the native FV1 integration-point velocity with the element-local
-- mixed solution. The steady flow equations carry no time mass and are
-- re-enforced at every backward-Euler transport level (quasi-steady flow).
transportSpace = approxSpace
-- Keep the pressure-stabilized FV1 velocity values while treating the
-- converged flow as a one-way coefficient in the transport Jacobian.  This is
-- essential on refined grids: concentration must not drive a flow correction.
stokesEq:set_velocity_ip_explicit(true)
transportStokesVelocity = stokesEq:velocity_ip()
outputStokesVelocity = stokesEq:velocity()

transportMembraneVelocity = DarcyVelocityLinker()
transportMembraneVelocity:set_permeability(membranePermeability * pressureScale)
transportMembraneVelocity:set_viscosity(kinematicViscosity)
transportMembraneVelocity:set_density(1.0)
transportMembraneVelocity:set_gravity(ConstUserVector(0.0))
transportMembraneVelocity:set_pressure_gradient(darcyMembraneEq:gradient())
transportMembraneVelocity:set_derivative_mask(1)

transportInterVelocity = DarcyVelocityLinker()
transportInterVelocity:set_permeability(interPermeability * pressureScale)
transportInterVelocity:set_viscosity(kinematicViscosity)
transportInterVelocity:set_density(1.0)
transportInterVelocity:set_gravity(ConstUserVector(0.0))
transportInterVelocity:set_pressure_gradient(darcyInterEq:gradient())
transportInterVelocity:set_derivative_mask(1)

transportDisc = DomainDiscretization(transportSpace)
transportDisc:add(stokesEq)
transportDisc:add(stokesInlet)
transportDisc:add(stokesOpenBoundary)
transportDisc:add(stokesRim)
transportDisc:add(stokesBJS)
transportDisc:add(darcyMembraneEq)
transportDisc:add(darcyInterEq)
transportDisc:add(pressureTrace)
transportDisc:add(darcyOuterPressure)
transportDisc:add(inactiveDarcy)
transportDisc:add(inactiveStokes)

for _, component in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    local coupledInterfaceTraction = NeumannBoundary(component[1], "fv1")
    coupledInterfaceTraction:add(
        darcyPressure * component[2], "Apical", "Lumen")
    transportDisc:add(coupledInterfaceTraction)

    local coupledOutletTraction = NeumannBoundary(component[1], "fv1")
    coupledOutletTraction:add(
        outletKinematicPressure * component[2], "Outlet", "Lumen")
    transportDisc:add(coupledOutletTraction)
end

transportEquations = {}
for _, soluteKey in ipairs(soluteOrder) do
    local parameter = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    local equations = {}

    equations.lumen = ConvectionDiffusionFV1(fields.lumen, "Lumen")
    equations.lumen:set_diffusion(parameter.diffusion)
    equations.lumen:set_velocity(transportStokesVelocity)
    equations.lumen:set_upwind(FullUpwind())
    transportDisc:add(equations.lumen)

    equations.membrane = ConvectionDiffusionFV1(
        fields.membrane, "Membrane")
    equations.membrane:set_diffusion(parameter.diffusion)
    equations.membrane:set_velocity(transportMembraneVelocity)
    equations.membrane:set_upwind(FullUpwind())
    transportDisc:add(equations.membrane)

    equations.inter = ConvectionDiffusionFV1(fields.inter, "Inter")
    equations.inter:set_diffusion(parameter.diffusion)
    equations.inter:set_velocity(transportInterVelocity)
    equations.inter:set_upwind(FullUpwind())
    transportDisc:add(equations.inter)

    if parameter.apicalPermeability > 0.0 then
        local leak = Leak({fields.lumen, fields.membrane})
        leak:set_permeability(parameter.apicalPermeability)
        local transport = MembraneTransportFV1("Apical", leak)
        transport:set_density_function(membraneTransportDensity)
        transportDisc:add(transport)
    end

    if parameter.basolateralPermeability > 0.0 then
        local leak = Leak({fields.membrane, fields.inter})
        leak:set_permeability(parameter.basolateralPermeability)
        local transport = MembraneTransportFV1("Basolateral", leak)
        transport:set_density_function(membraneTransportDensity)
        transportDisc:add(transport)
    end

    -- FV1 convection needs the open-boundary c*u.n term explicitly.
    local outletFlux = ScaleAddLinkerVector()
    outletFlux:add(equations.lumen:value(), transportStokesVelocity)
    local outletBoundary = NeumannBoundary(fields.lumen, "fv1")
    outletBoundary:add(outletFlux, "Outlet", "Lumen")
    transportDisc:add(outletBoundary)

    transportEquations[soluteKey] = equations
end

concentrationBoundary = DirichletBoundary()
for _, soluteKey in ipairs(soluteOrder) do
    local parameter = soluteParameters[soluteKey]
    concentrationBoundary:add(
        parameter.inlet, transportFields[soluteKey].lumen, "Inlet")
end
transportDisc:add(concentrationBoundary)
transportDisc:add(OneSideP1Constraints())

transportSolution = solution:clone()
for _, soluteKey in ipairs(soluteOrder) do
    local parameter = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    Interpolate(parameter.inlet, transportSolution, fields.lumen)
    Interpolate(parameter.cellInitial, transportSolution, fields.membrane)
    Interpolate(parameter.interInitial, transportSolution, fields.inter)
end
transportDisc:adjust_solution(transportSolution)

initialMasses = {}
for _, soluteKey in ipairs(soluteOrder) do
    local fields = transportFields[soluteKey]
    initialMasses[soluteKey] =
        Integral(transportSolution, fields.lumen, "Lumen") +
        Integral(transportSolution, fields.membrane, "Membrane") +
        Integral(transportSolution, fields.inter, "Inter")
end

--------------------------------------------------------------------------------
-- Scalable nonlinear solver
--------------------------------------------------------------------------------

-- Full-system multigrid block smoother.  Uzawa isolates the Stokes pressure;
-- ILU acts on the remaining flow and four algebraically independent
-- three-compartment transport blocks while preserving their matrix sparsity.
transportRegularInverse = ILU()
transportRegularInverse:set_damp(0.8)
transportRegularInverse:set_inversion_eps(1.0e-30)
transportPressureInverse = ILU()
transportPressureInverse:set_damp(1.0)
transportPressureInverse:set_inversion_eps(1.0e-30)
transportUzawa = UzawaBase({"p_stokes"})
transportUzawa:set_forward_iter(transportRegularInverse)
transportUzawa:set_schur_iter(transportPressureInverse)
transportUzawa:set_backward_iter(transportRegularInverse)
transportUzawa:set_schur_operator_update(
    schurMassOp, 1.0 / kinematicViscosity)

util.solver.defaults.approxSpace = transportSpace
transportSolverConfig = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 30,
        absolute = 1.0e-10, reduction = 1.0e-8, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg", approxSpace = transportSpace,
            discretization = transportDisc, baseLevel = 0,
            smoother = transportUzawa,
            cycle = "V", preSmooth = 2, postSmooth = 2,
            -- Rediscretize the complete coupled operator on each level.  RAP
            -- mixes the differently supported 17-field blocks and destroys
            -- the one-way block structure on the coarse level.
            rap = false, transfer = "std", baseSolver = "superlu",
            gatheredBaseSolverIfAmbiguous = true
        },
        convCheck = {
            type = "standard", iterations = 400,
            absolute = 1.0e-12, reduction = 1.0e-9, verbose = true
        }
    }
}
transportSolver = util.solver.CreateSolver(transportSolverConfig)

--------------------------------------------------------------------------------
-- Compartment time-series output
--------------------------------------------------------------------------------

outputTimes = {}
lastOutputStep = 0
compartmentOutputs = {
    {name = "lumen", subset = "Lumen", velocity = transportStokesVelocity,
     writer = VTKOutput()},
    {name = "membrane", subset = "Membrane",
     velocity = transportMembraneVelocity, writer = VTKOutput()},
    {name = "inter", subset = "Inter", velocity = transportInterVelocity,
     writer = VTKOutput()},
}
for _, compartment in ipairs(compartmentOutputs) do
    compartment.writer:clear_selection()
    for _, soluteKey in ipairs(soluteOrder) do
        local parameter = soluteParameters[soluteKey]
        compartment.writer:select(
            transportFields[soluteKey][compartment.name], parameter.label)
    end
    compartment.writer:select_element(
        compartment.velocity, "velocity_mm_s")
end

function WriteSerialPVTU(prefix, step, time)
    if ProcRank() ~= 0 or NumProcs() ~= 1 then return end
    local fileName = string.format("%s_t%04d.pvtu", prefix, step)
    local file = assert(io.open(fileName, "w"),
        "Cannot create PVTU file: " .. fileName)
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n')
    file:write('<VTKFile type="PUnstructuredGrid" version="0.1" ',
               'byte_order="LittleEndian">\n')
    file:write(string.format('  <Time timestep="%.17g"/>\n', time))
    file:write('  <PUnstructuredGrid GhostLevel="0">\n')
    file:write('    <PPoints><PDataArray type="Float32" ',
               'NumberOfComponents="3"/></PPoints>\n')
    file:write('    <PPointData>\n')
    for _, soluteKey in ipairs(soluteOrder) do
        local parameter = soluteParameters[soluteKey]
        file:write(string.format(
            '      <PDataArray type="Float32" Name="%s" ', parameter.label),
            'NumberOfComponents="1"/>\n')
    end
    file:write('    </PPointData>\n')
    file:write('    <PCellData><PDataArray type="Float32" ',
               'Name="velocity_mm_s" NumberOfComponents="3"/></PCellData>\n')
    file:write(string.format('    <Piece Source="%s_t%04d.vtu"/>\n',
                             basename, step))
    file:write('  </PUnstructuredGrid>\n</VTKFile>\n')
    file:close()
end

function WriteCompartmentPVD(prefix)
    if ProcRank() ~= 0 then return end
    local file = assert(io.open(prefix .. ".pvd", "w"),
        "Cannot create PVD file: " .. prefix .. ".pvd")
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n')
    file:write('<VTKFile type="Collection" version="0.1" ',
               'byte_order="LittleEndian">\n  <Collection>\n')
    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            file:write(string.format(
                '    <DataSet timestep="%.17g" group="" part="0" ' ..
                'file="%s_t%04d.pvtu"/>\n',
                outputTimes[step], basename, step))
        end
    end
    file:write('  </Collection>\n</VTKFile>\n')
    file:close()
end

function WriteTransportOutput(currentSolution, step, time)
    for _, compartment in ipairs(compartmentOutputs) do
        local prefix = outputName .. "_" .. compartment.name
        compartment.writer:print_subsets(
            prefix, currentSolution, compartment.subset, step, time)
        WriteSerialPVTU(prefix, step, time)
    end
    outputTimes[step] = time
    lastOutputStep = math.max(lastOutputStep, step)
    for _, compartment in ipairs(compartmentOutputs) do
        WriteCompartmentPVD(outputName .. "_" .. compartment.name)
    end
end

WriteTransportOutput(transportSolution, 0, 0.0)

timeDisc = ThetaTimeStep(transportDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(transportSolver)
timeIntegrator:set_time_step(dt)
transportSolver:init(AssembledOperator(transportDisc))
transportSolver:prepare(transportSolution)

transportObserver = LuaCallbackObserver()
function StandaloneTransportOutputCallback(step, time, currentDt)
    if step % outputEverySteps == 0 or time >= endTime - 1.0e-12 then
        WriteTransportOutput(
            transportObserver:get_current_solution(), step, time)
    end
    return 1
end
transportObserver:set_callback("StandaloneTransportOutputCallback")
timeIntegrator:attach_observer(transportObserver)

print("Starting transient transport solve...")
timeIntegrator:apply(
    transportSolution, endTime, transportSolution, 0.0)

for _, soluteKey in ipairs(soluteOrder) do
    local parameter = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    local finalMass =
        Integral(transportSolution, fields.lumen, "Lumen") +
        Integral(transportSolution, fields.membrane, "Membrane") +
        Integral(transportSolution, fields.inter, "Inter")
    print(string.format(
        "%s total mass integral: initial=%.17e final=%.17e change=%.17e",
        parameter.label, initialMasses[soluteKey], finalMass,
        finalMass - initialMasses[soluteKey]))
end

print("Transport output prefix: " .. outputName)
print("STANDALONE_FLOW_TRANSPORT: COMPLETE")

