------------------------------------------------------------------------------------------
--  Steady Stokes--Darcy flow in the full fourth nephron trace, FE Darcy
--
--  Author:   Baleeg Alhilali
--  Date:     2026-11-26 (restyled; physics and solver from the 2026-11-24 version)
--  Based on: code/flow/nephron_stokes_darcy_clean.lua
--
--  Lumen:        Stokes FV1
--  Membrane:     Darcy FE
--  Interstitium: Darcy FE
--  Interface:    pressure trace, Darcy normal traction and BJS slip on Apical
--
--  Run (8 ranks, Model root):
--  mpirun -np 8 <UG4>/bin/ugshell -ex code/flow/nephron_flow_fourth_pad100.lua -numRefs 2
--  Results: docs/FOURTH_TRACE_FULL_FLOW_2026-11-24.md
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------
-- Initiation
------------------------------------------------------------------------------------------

totalStartTime = GetClockS()

-- Model root from this script's location (Model/code/flow/<script>.lua)
scriptPath = debug.getinfo(1, "S").source
if scriptPath:sub(1, 1) == "@" then scriptPath = scriptPath:sub(2) end
if scriptPath:sub(1, 1) ~= "/" then
    scriptPath = CurrentWorkingDirectory() .. "/" .. scriptPath
end
MODEL_ROOT = scriptPath:match("^(.*)/code/flow/[^/]+$")
if MODEL_ROOT == nil then
    print("Expected this script inside Model/code/flow: " .. scriptPath)
    exit()
end

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

dim = 3
InitUG(dim, AlgebraType("CPU", 1))

------------------------------------------------------------------------------------------
-- Parameters
------------------------------------------------------------------------------------------

-- Units: mm, s; velocity mm/s, kinematic pressure mm^2/s^2, permeability mm^2

-- Fluid parameters
kinematicViscosity = 1.0                         -- [mm^2/s]

-- Permeabilities
membranePermeability = 7.0e-12                   -- [mm^2]  ~9 % reabsorption
interPermeability    = 1.0e-11                   -- [mm^2]
bjsPermeability      = 1.0e-13                   -- [mm^2]  Beavers-Joseph-Saffman

-- Interface parameters
bjsAlpha        = 0.2                            -- [-]
normalImpedance = 1.0e3                          -- [mm/s]

-- Boundary pressures
outerKinematicPressure  = 6.132812e5             -- [mm^2/s^2]  OuterWall (613 Pa)
outletKinematicPressure = 2.6664474e6            -- [mm^2/s^2]  Outlet

-- Inflow
inletFlowRate = 5.0e-4                           -- [mm^3/s]  30 nL/min

-- Numerical parameters
pressureScale                         = 1.0e4    -- [-]  p_darcy_scaled = p / pressureScale
firstDistributionLevel                = 0        -- [-]
interfaceConservationTolerancePercent = 1.0e-1   -- [%]  output gate

-- Mesh-dependent parameters (fourth_trace_full_pad100.ugx, inlet on the x-min wall)
inletCenter    = {0.360781065, 0.393057145, 0.162855950}  -- [mm]
inletDirection = {1.0, 0.0, 0.0}                          -- [-]
inletRadius    = 1.0e-2                                   -- [mm]

-- Poiseuille inlet: u_mean = Q / (pi R^2), u_max = 2 u_mean
meanInletVelocity    = inletFlowRate / (math.pi * inletRadius * inletRadius)  -- [mm/s]
maximumInletVelocity = 2.0 * meanInletVelocity                                -- [mm/s]

------------------------------------------------------------------------------------------
-- Terminal-dependent parameters
------------------------------------------------------------------------------------------

gridName = util.GetParam("-grid",
    MODEL_ROOT .. "/inputs/fourth_trace_full_pad100/fourth_trace_full_pad100.ugx",
    "grid file")
numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")

gridBaseName = gridName:match("([^/\\]+)%.%w+$")
outputName = util.GetParam("-output",
    MODEL_ROOT .. "/runs/flow_fe_" .. gridBaseName .. "_" .. os.date("%d%b_%Hh"),
    "VTK output prefix")

------------------------------------------------------------------------------------------
-- Functions
------------------------------------------------------------------------------------------

-- Parabolic (Poiseuille) inlet velocity, zero on and outside the inlet rim
function InletVelocity(x, y, z, t)
    local c, d = inletCenter, inletDirection
    local dx, dy, dz = x - c[1], y - c[2], z - c[3]
    local axial = dx * d[1] + dy * d[2] + dz * d[3]
    local rx, ry, rz = dx - axial * d[1], dy - axial * d[2], dz - axial * d[3]
    local radiusSquared = rx * rx + ry * ry + rz * rz
    local profile = math.max(0.0, 1.0 - radiusSquared / (inletRadius * inletRadius))
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return maximumInletVelocity * profile * d[1],
        maximumInletVelocity * profile * d[2],
        maximumInletVelocity * profile * d[3]
end

-- Element-wise x, y, z components of a vector field (for VTK output)
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
    output:select_nodal("p_darcy_scaled", "pressure_kinematic")
    output:select_element(velocity, "velocity_mm_s")
    output:select_element(components[1], "u")
    output:select_element(components[2], "v")
    output:select_element(components[3], "w")
    output:print_subsets(outputName .. suffix, physicalSolution, subset, 0, 0.0)
end

------------------------------------------------------------------------------------------
-- Discretization
------------------------------------------------------------------------------------------

print("Steady nephron Stokes--Darcy flow (FE Darcy)")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())
print(string.format(
    "Inlet: Q=%.6e mm^3/s (%.6g nL/min), R=%.6e mm, u_mean=%.6e mm/s, u_max=%.6e mm/s",
    inletFlowRate, inletFlowRate * 6.0e4, inletRadius,
    meanInletVelocity, maximumInletVelocity))

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)

balancer.firstDistLvl = firstDistributionLevel
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, numRefs)

print(dom:domain_info():to_string())

allSubsets = "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall," ..
    "MembraneWall,Inlet,Outlet,TripleJunction"
velocitySubsets = {"Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"}
stokesPressureSubsets = {"Lumen", "Membrane", "Apical", "Basolateral",
    "MembraneWall", "Inlet", "Outlet", "TripleJunction"}

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct({"u", "v", "w"}, "Lagrange", 1, velocitySubsets)
approxSpace:add_fct({"p_stokes"}, "Lagrange", 1, stokesPressureSubsets)
approxSpace:add_fct("p_darcy_scaled", "Lagrange", 1, allSubsets)
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic("component")

--------------------------------------------------------------------------------
-- Stokes equation (Lumen)
--------------------------------------------------------------------------------

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
stokesEq:set_stokes(true)
stokesEq:set_laplace(false)
stokesEq:set_kinematic_viscosity(kinematicViscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_upwind("lps")
stokesEq:set_peclet_blend(false)
stokesEq:set_stabilization("fields", "cor")
stokesEq:set_pac_upwind(false)

-- Inlet: parabolic velocity profile
stokesInlet = NavierStokesInflow(stokesEq)
stokesInlet:add("InletVelocity", "Inlet")

-- Outlet: no normal viscous stress; the outlet pressure enters as traction below
stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
stokesOpenBoundary:add("Outlet")

-- Lumen rim: no slip
stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

-- Apical: BJS tangential friction and normal impedance
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsAlpha * kinematicViscosity / math.sqrt(bjsPermeability))
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalImpedance)

-- Stokes pressure is inactive outside the lumen
inactiveStokes = DirichletBoundary()
inactiveStokes:add(0.0, "p_stokes", "Membrane,Basolateral,MembraneWall")

--------------------------------------------------------------------------------
-- Darcy equation (Membrane and Interstitium)
--------------------------------------------------------------------------------

darcyMembraneEq = ConvectionDiffusionFE("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(membranePermeability / kinematicViscosity * pressureScale)

darcyInterEq = ConvectionDiffusionFE("p_darcy_scaled", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(interPermeability / kinematicViscosity * pressureScale)

-- OuterWall: fixed interstitial pressure
darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(outerKinematicPressure / pressureScale, "p_darcy_scaled", "OuterWall")

-- Darcy pressure is inactive inside the lumen
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet")

--------------------------------------------------------------------------------
-- Interface coupling (Apical)
--------------------------------------------------------------------------------

-- Pressure continuity p_stokes = pressureScale * p_darcy_scaled; merges the
-- Stokes and Darcy flux rows so the normal flux is conserved
pressureTrace = PressureTraceConstraint("p_stokes", "p_darcy_scaled", "Apical,TripleJunction")
pressureTrace:set_pressure_scale(pressureScale)

--------------------------------------------------------------------------------
-- Global discretization
--------------------------------------------------------------------------------

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

-- Normal traction: Darcy pressure on Apical, outlet pressure on Outlet
darcyPressure = pressureScale * darcyMembraneEq:value()
for _, component in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    local interfaceTraction = NeumannBoundary(component[1], "fv1")
    interfaceTraction:add(darcyPressure * component[2], "Apical", "Lumen")
    domainDisc:add(interfaceTraction)

    local outletTraction = NeumannBoundary(component[1], "fv1")
    outletTraction:add(outletKinematicPressure * component[2], "Outlet", "Lumen")
    domainDisc:add(outletTraction)
end

------------------------------------------------------------------------------------------
-- Solver
------------------------------------------------------------------------------------------

-- Uzawa smoother: line Gauss-Seidel on velocity, ILU on the Stokes pressure Schur
-- complement approximated by the pressure mass matrix / viscosity
lineGS = LineGaussSeidel(approxSpace)
lineGS:set_num_steps(1, 1, 1, 1, 1, 1)

pressureILU = ILU()
pressureILU:set_damp(1.0)
pressureILU:set_inversion_eps(1.0e-30)

schurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
schurMassEq:set_mass_scale(1.0)
schurMassEq:set_diffusion(0.0)
schurMassDisc = DomainDiscretization(approxSpace)
schurMassDisc:add(schurMassEq)
schurMassDisc:add(OneSideP1Constraints())
schurMassOp = AssembledLinearOperator(schurMassDisc)

uzawa = UzawaBase({"p_stokes"})
uzawa:set_forward_iter(lineGS)
uzawa:set_schur_iter(pressureILU)
uzawa:set_backward_iter(lineGS)
uzawa:set_schur_operator_update(schurMassOp, 1.0 / kinematicViscosity)

solverDesc = {
    type = "bicgstab",                  -- ["bicgstab", "gmres", "cg", "linear"]
    precond = {
        type       = "gmg",             -- ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
        approxSpace    = approxSpace,
        discretization = domainDisc,
        rap        = true,              -- Galerkin coarse operators
        smoother   = uzawa,             -- ["ilu", "ilut", "jac", "gs", "sgs"] or a Lua smoother object
        cycle      = "V",               -- ["V", "F", "W"]
        preSmooth  = 2,
        postSmooth = 2,
        baseSolver = "superlu",         -- ["superlu", "lu", "ilu", ...]; coarse level 0 only, to be replaced
        baseLevel  = 0,
        gatheredBaseSolverIfAmbiguous = true,
        transfer = {
            type = "std",
            restrictionDamp  = 1.0,
            prolongationDamp = 1.0,
            enableP1LagrangeOptimization = true
        }
    },
    convCheck = {
        type       = "standard",
        iterations = 400,               -- maximum number of iterations
        absolute   = 1.0e-14,           -- absolute defect to be reached
        reduction  = 1.0e-8,            -- defect reduction to be reached
        verbose    = true
    }
}

util.solver.defaults.approxSpace = approxSpace
solver = util.solver.CreateLinearSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Initial conditions
------------------------------------------------------------------------------------------

solution = GridFunction(approxSpace)
solution:set(0.0)
domainDisc:adjust_solution(solution)
solution:set_consistent_storage_type()

------------------------------------------------------------------------------------------
-- Solve (steady)
------------------------------------------------------------------------------------------

matrix = AssembledLinearOperator(domainDisc)
rightHandSide = solution:clone()
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

------------------------------------------------------------------------------------------
-- Output
------------------------------------------------------------------------------------------

-- Assembled interface rows, captured before the pressure trace merges them:
-- the exact algebraic normal-flux conservation check
interfaceDefect = solution:clone()
pressureTrace:set_interface_defect_dump(outputName .. "_interface")
pressureTrace:capture_next_interface_defect()
domainDisc:assemble_defect(interfaceDefect, solution)
interfaceDefectSums = pressureTrace:interface_defect_sums()
stokesInterfaceResidualFlux = interfaceDefectSums[1]
darcyInterfaceResidualFlux = interfaceDefectSums[2]
interfaceFluxMismatch = interfaceDefectSums[3]

physicalSolution = solution:clone()
ScaleGF(physicalSolution, solution, {1.0, 1.0, 1.0, 1.0, pressureScale})
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

inletFlow = IntegralNormalComponentOnManifold(stokesVelocity, physicalSolution, "Inlet", "Lumen")
outletFlow = IntegralNormalComponentOnManifold(stokesVelocity, physicalSolution, "Outlet", "Lumen")
apicalFlow = IntegralNormalComponentOnManifold(stokesVelocity, physicalSolution, "Apical", "Lumen")
darcyApicalFlow = IntegralNormalComponentOnManifold(
    darcyMembraneVelocity, physicalSolution, "Apical", "Membrane")

inletPressureMean = Integral(physicalSolution, "p_stokes", "Inlet")
    / Integral(1.0, physicalSolution, "Inlet")
outletPressureMean = Integral(physicalSolution, "p_stokes", "Outlet")
    / Integral(1.0, physicalSolution, "Outlet")
pressureDrop = inletPressureMean - outletPressureMean

interfaceFluxScale = math.max(math.abs(stokesInterfaceResidualFlux),
    math.abs(darcyInterfaceResidualFlux), 1.0e-30)
interfaceFluxMismatchPercent = 100.0 * math.abs(interfaceFluxMismatch) / interfaceFluxScale

flowBalance = inletFlow + outletFlow + apicalFlow
flowBalancePercent = 100.0 * math.abs(flowBalance) / math.abs(inletFlow)

print(string.format("Flow [mm^3/s]: inlet=%.6e, outlet=%.6e, apical=%.6e",
    inletFlow, outletFlow, apicalFlow))
print(string.format("Mean pressure [mm^2/s^2]: inlet=%.12e, outlet=%.12e, drop=%.12e",
    inletPressureMean, outletPressureMean, pressureDrop))
print(string.format("Lumen balance: %.6e mm^3/s (%.4f%%)", flowBalance, flowBalancePercent))
print(string.format("Assembled interface rows [mm^3/s]: Stokes=%.6e, Darcy=%.6e",
    stokesInterfaceResidualFlux, darcyInterfaceResidualFlux))
print(string.format("Assembled interface conservation: sum=%.6e mm^3/s (%.6g%%)",
    interfaceFluxMismatch, interfaceFluxMismatchPercent))
print(string.format("Reconstructed Apical fluxes [mm^3/s]: Stokes=%.6e, Darcy=%.6e",
    apicalFlow, darcyApicalFlow))

if interfaceFluxMismatchPercent > interfaceConservationTolerancePercent then
    print(string.format("Interface conservation failed: %.6g%% > %.6g%%",
        interfaceFluxMismatchPercent, interfaceConservationTolerancePercent))
    exit()
end

lumenOutput = VTKOutput()
lumenOutput:clear_selection()
lumenOutput:select_nodal({"u", "v", "w"}, "velocity_mm_s")
lumenOutput:select_nodal("p_stokes", "pressure_kinematic")
lumenOutput:print_subsets(outputName .. "_lumen", physicalSolution, "Lumen", 0, 0.0)

WriteDarcyOutput("_membrane", "Membrane", darcyMembraneVelocity)
WriteDarcyOutput("_inter", "Inter", darcyInterVelocity)

print(string.format("Solver time: %.3f s", solveTime))
print(string.format("Total time:  %.3f s", GetClockS() - totalStartTime))
print("Done")

------------------------------------------------------------------------------------------
-- End of File
------------------------------------------------------------------------------------------
