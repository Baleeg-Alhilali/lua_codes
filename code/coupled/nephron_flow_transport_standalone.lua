----------------------------------------------------------
--
--  Standalone quasi-steady nephron flow and multispecies leak transport
--
--  Flow (solved once, steady):
--    Lumen:        Stokes FV1
--    Membrane:     Darcy FE (or FV1, see darcyScheme)
--    Interstitium: Darcy FE (or FV1)
--    Interface:    pressure trace, Darcy traction, BJS friction;
--                  NoNormalStressOutflow on Outlet only (2026-11-24 fix)
--
--  Transport (transient, backward Euler):
--    Na+, K+, Cl-, urea advection--diffusion in Lumen, Membrane, Inter
--    Apical and Basolateral: passive leaks (GHK for Na and Cl at fixed
--    potentials), basolateral Na/K pump, solvent drag (1 - sigma) J_v c_up
--    with the conservative membrane water flux J_v
--    Reabsorption bookkeeping [pmol] printed and written to
--    <output>_reabsorption.csv
--
--  Coupling: flow and concentrations share one approximation space so that
--  transport convects with the same FV1 face fluxes the flow conserved
--  (Stokes velocity_ip, Darcy pressure gradient).  Block solve per time step
--  (block Gauss-Seidel; exact for the one-way coupling): (1) the
--  Stokes--Darcy block with its own GMG, concentrations held; (2) the
--  transport block with its own GMG, flow DoFs held at the new flow.
--
--  This is a software/transport test stage, not a physiological model.
--  See docs/LEAK_TRANSPORT_BLOCK_SOLVE_2026-11-26.md.
--
--  Run from any directory, e.g. (8 ranks):
--    mpirun -np 8 <UG4>/bin/ugshell -ex code/coupled/nephron_flow_transport_standalone.lua \
--        -numRefs 0 -endTime 20 -dt 0.1 -outputEverySteps 10 -output <prefix>
--
----------------------------------------------------------

-- Resolve the Model root from this script's location.
local __script = debug.getinfo(1, "S").source
if __script:sub(1, 1) == "@" then __script = __script:sub(2) end
if __script:sub(1, 1) ~= "/" then
    __script = CurrentWorkingDirectory() .. "/" .. __script
end
MODEL_ROOT = assert(__script:match("^(.*)/code/coupled/[^/]+$"),
    "Expected this script inside Model/code/coupled: " .. __script)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = 3

numRefs = util.GetParamNumber("-numRefs", 0, "number of refinements")
gridName = util.GetParam("-grid", MODEL_ROOT ..
    "/mesh_generation/meshes/2026-11-26_fourth_best_mesh/fourth_trace_best_mesh.ugx",
    "coarse projector-backed UGX (meshes live in Model/mesh_generation/meshes)")
outputName = util.GetParam(
    "-output", MODEL_ROOT .. "/Results/" .. os.date("%Y-%m-%d") .. "_leak_transport/leak_transport_r" .. numRefs,
    "VTK output prefix")

-- Transient transport run control
endTime = util.GetParamNumber("-endTime", 30.0)
dt = util.GetParamNumber("-dt", 0.1)
outputEverySteps = util.GetParamNumber("-outputEverySteps", 10)

-- Internal unit system: millimetres and seconds.
-- Velocity: mm/s; kinematic pressure: mm^2/s^2; permeability: mm^2;
-- concentration: mM; membrane permeability: mm/s.

-- Inlet (fourth_trace_best_mesh.ugx: 10 um hex Inter shell, x-min wall).  Poiseuille profile with
-- u_mean = Q/(pi R^2), u_max = 2 u_mean.
inletCenterX = 0.338568
inletCenterY = 0.400919
inletCenterZ = 0.163024
inletDirectionX = 1.0
inletDirectionY = 0.0
inletDirectionZ = 0.0
inletRadius = 1.0e-2
inletFlowRate = 5.0e-4                    -- [mm^3/s] = 30 nL/min
meanInletVelocity = inletFlowRate / (math.pi * inletRadius * inletRadius)
maximumInletVelocity = 2.0 * meanInletVelocity

-- Fluid and porous-medium parameters (docs/FOURTH_TRACE_FULL_FLOW_2026-11-24.md)
kinematicViscosity = 1.0                  -- [mm^2/s]
membranePermeability = 7.0e-12            -- [mm^2] (~8.8% reabsorption)
interPermeability = 1.0e-11               -- [mm^2]
bjsPermeability = 1.0e-13                 -- [mm^2]
bjsAlpha = 0.2
normalImpedance = 1.0e3                   -- [mm/s]
pressureScale = 1.0e4
outerKinematicPressure = 6.132812e5       -- [mm^2/s^2] OuterWall
outletKinematicPressure = 2.6664474e6     -- [mm^2/s^2] Outlet

-- Darcy discretization: "FE" (validated BiCGStab + Uzawa/line-GS GMG flow
-- solver) or "FV1" (box-conservative, GMRES(100) + Uzawa/ILU GMG).
darcyScheme = "FE"

-- Parallel distribution and validation
firstDistributionLevel = 0
interfaceConservationTolerancePercent = 1.0e-1
massBalanceTolerancePercent = 1.0e-6      -- per-step discrete balance gate

-- Solutes.  Free-solution diffusion in all compartments (as the validated
-- urea model).  Leak permeabilities [mm/s] act on the smooth face area, so no
-- radius factor is applied; a zero permeability means no leak.  TEST STAGE
-- (2026-11-26): Na leaks only apically and leaves the cell by the
-- basolateral Na/K pump; K leaks only basolaterally (plus the pump); Cl and
-- urea leak across both membranes.  Weinstein values where nonzero; the
-- 2026-11-22 project assumptions are marked "project".
-- Reflection coefficients sigma [0,1]: the solvent-drag flux across a
-- membrane is (1 - sigma) J_v c_upstream (sigma = 1: no drag, the solute is
-- left behind when water is reabsorbed; sigma = 0: it moves with the water).
soluteParameters = {
    na = {
        label = "Na", diffusion = 1.33e-3, inlet = 140.0,
        valency = 1, ghk = false,            -- [-] GHK leaks
        cellInitial = 19.6, interInitial = 0,
        outerWall = 0,                    -- [mM] OuterWall bath
        apicalPermeability = 1.00e-3,       -- project
        basolateralPermeability = 0.0,      -- pump pathway only
        reflectionApical = 0.0,             -- [-] solvent drag, test value
        reflectionBasolateral = 0.0,        -- [-] solvent drag, test value
    },
    k = {
        label = "K", diffusion = 1.96e-3, inlet = 4.9,
        valency = 1, ghk = false,           -- [-] linear leak
        cellInitial = 2.0, interInitial = 4.66,
        outerWall = 4.66,                   -- [mM] OuterWall bath
        apicalPermeability = 0.0,           -- no apical K leak
        basolateralPermeability = 7.40e-4,
        reflectionApical = 0.0,             -- [-] solvent drag, test value
        reflectionBasolateral = 0.0,        -- [-] solvent drag, test value
    },
    cl = {
        label = "Cl", diffusion = 2.03e-3, inlet = 113.2,
        valency = -1, ghk = false,           -- [-] GHK leaks
        cellInitial = 50.0, interInitial = 0,
        outerWall = 0.0,                    -- [mM] OuterWall bath
        apicalPermeability = 8.00e-4,       -- project
        basolateralPermeability = 8.00e-4,  -- project
        reflectionApical = 0.0,             -- [-] solvent drag, test value
        reflectionBasolateral = 0.0,        -- [-] solvent drag, test value
    },
    urea = {
        label = "Urea", diffusion = 1.67e-3, inlet = 5.0,
        valency = 0, ghk = false,           -- [-] linear leak
        cellInitial = 1.0, interInitial = 0,
        outerWall = 0.0,                    -- [mM] OuterWall bath
        apicalPermeability = 3.78e-4,
        basolateralPermeability = 3.70e-4,
        reflectionApical = 0.0,             -- [-] solvent drag, test value
        reflectionBasolateral = 0.0,        -- [-] solvent drag, test value
    },
}
soluteOrder = {"na", "k", "cl", "urea"}

-- Basolateral Na/K-ATPase (neuro_collection NaKPump), 3 Na out : 2 K in,
-- J_cycle = Jmax [Na_i/(Na_i+K_Na)]^3 [K_o/(K_o+K_K)]^2 (2026-11-22 note:
-- 10.8 nmol/(s cm^2) Na = 0.108 mM mm/s, i.e. Jmax = 0.036 mM mm/s).
pumpEnabled = true
pumpMaxCycleFlux = 0.036                  -- [mM mm/s]
pumpKNa = 3.516                           -- [mM]
pumpKK = 0.8584                           -- [mM]

-- Electrical potentials for the GHK leaks (fixed, no potential equation).
-- Cell (Membrane compartment) at -70 mV relative to Lumen and Inter.
faradayConstant = 96485.0                 -- [C/mol]
gasConstant = 8.31451                     -- [J/(mol K)]
temperature = 310.15                      -- [K]
compartmentPotential = {lumen = 0.0, membrane = -0.050, inter = 0.0}  -- [V]
compartmentOrder = {"lumen", "membrane", "inter"}

assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0 and dt <= endTime, "-dt must be in (0, endTime]")
assert(outputEverySteps >= 1 and
       outputEverySteps == math.floor(outputEverySteps),
       "-outputEverySteps must be a positive integer")

--------------------------------------------------------------------------------
-- Initialization
--------------------------------------------------------------------------------

local function EnsureOutputDirectory(prefix)
    local directory = string.match(prefix, "^(.*)/[^/]*$")
    if directory == nil or directory == "" then return end
    local current = string.sub(directory, 1, 1) == "/" and "/" or ""
    for part in string.gmatch(directory, "[^/]+") do
        if current == "/" or current == "" then
            current = current .. part
        else
            current = current .. "/" .. part
        end
        if not DirectoryExists(current) then
            assert(CreateDirectory(current) or DirectoryExists(current),
                   "Could not create output directory: " .. current)
        end
    end
end

totalStartTime = GetClockS()
InitUG(dim, AlgebraType("CPU", 1))

if ProcRank() == 0 then EnsureOutputDirectory(outputName) end
SynchronizeProcesses()

print("Standalone quasi-steady nephron flow (" .. darcyScheme .. " Darcy) and leak transport")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())
print(string.format(
    "Inlet: Q=%.6e mm^3/s, R=%.6e mm, u_mean=%.6e mm/s",
    inletFlowRate, inletRadius, meanInletVelocity))
print(string.format("Transport: 0 to %.6g s, dt=%.6g s, output every %d steps",
                    endTime, dt, outputEverySteps))
for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    print(string.format(
        "%s: inlet=%.6g cell0=%.6g inter0=%.6g outerWall=%.6g mM; " ..
        "D=%.6g mm^2/s; Papical=%.6g Pbaso=%.6g mm/s",
        p.label, p.inlet, p.cellInitial, p.interInitial, p.outerWall, p.diffusion,
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
print(dom:domain_info():to_string())

allSubsets =
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall," ..
    "MembraneWall,Inlet,Outlet,TripleJunction"
velocitySubsets = "Lumen,Apical,Inlet,Outlet,TripleJunction"
stokesPressureSubsets =
    "Lumen,Membrane,Apical,Basolateral,MembraneWall,Inlet,Outlet,TripleJunction"

transportFieldSubsets = {
    lumen = "Lumen,Apical,Inlet,Outlet,TripleJunction",
    membrane = "Membrane,Apical,Basolateral,MembraneWall,TripleJunction",
    inter = "Inter,Basolateral,OuterWall",
}
compartmentVolume = {lumen = "Lumen", membrane = "Membrane", inter = "Inter"}

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("u,v,w", "Lagrange", 1, velocitySubsets)
approxSpace:add_fct("p_stokes", "Lagrange", 1, stokesPressureSubsets)
approxSpace:add_fct("p_darcy_scaled", "Lagrange", 1, allSubsets)
transportFields = {}
for _, soluteKey in ipairs(soluteOrder) do
    transportFields[soluteKey] = {}
    for _, compartment in ipairs(compartmentOrder) do
        local field = soluteKey .. "_" .. compartment
        transportFields[soluteKey][compartment] = field
        approxSpace:add_fct(field, "Lagrange", 1,
                            transportFieldSubsets[compartment])
    end
end
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Flow discretization
--------------------------------------------------------------------------------

function InletVelocity(x, y, z, t)
    local dx = x - inletCenterX
    local dy = y - inletCenterY
    local dz = z - inletCenterZ
    local axial = dx * inletDirectionX + dy * inletDirectionY
        + dz * inletDirectionZ
    local rx = dx - axial * inletDirectionX
    local ry = dy - axial * inletDirectionY
    local rz = dz - axial * inletDirectionZ
    local r2 = rx * rx + ry * ry + rz * rz
    local profile = math.max(0.0, 1.0 - r2 / (inletRadius * inletRadius))
    if r2 >= inletRadius * inletRadius * (1.0 - 1.0e-6) then profile = 0.0 end
    return maximumInletVelocity * profile * inletDirectionX,
        maximumInletVelocity * profile * inletDirectionY,
        maximumInletVelocity * profile * inletDirectionZ
end

scaledMembraneMobility = membranePermeability / kinematicViscosity * pressureScale
scaledInterMobility = interPermeability / kinematicViscosity * pressureScale
bjsFriction = bjsAlpha * kinematicViscosity / math.sqrt(bjsPermeability)

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
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

-- Outlet only: on Apical the continuity term of this condition was a mass
-- sink after the pressure-trace merge (2026-11-24).
stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
stokesOpenBoundary:add("Outlet")

stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalImpedance)

assert(darcyScheme == "FE" or darcyScheme == "FV1", "darcyScheme: FE or FV1")
local DarcyEquation = darcyScheme == "FE" and ConvectionDiffusionFE
                      or ConvectionDiffusionFV1
darcyMembraneEq = DarcyEquation("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneMobility)

darcyInterEq = DarcyEquation("p_darcy_scaled", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterMobility)

darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(outerKinematicPressure / pressureScale,
                       "p_darcy_scaled", "OuterWall")
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet")
inactiveStokes = DirichletBoundary()
inactiveStokes:add(0.0, "p_stokes", "Membrane,Basolateral,MembraneWall")

pressureTrace = PressureTraceConstraint(
    "p_stokes", "p_darcy_scaled", "Apical,TripleJunction")
pressureTrace:set_pressure_scale(pressureScale)

-- Concentrations are fixed (zero) during the steady flow solve only.
steadyTransportBoundary = DirichletBoundary()
for _, soluteKey in ipairs(soluteOrder) do
    for _, compartment in ipairs(compartmentOrder) do
        steadyTransportBoundary:add(0.0, transportFields[soluteKey][compartment],
                                    transportFieldSubsets[compartment])
    end
end

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
domainDisc:add(steadyTransportBoundary)
domainDisc:add(OneSideP1Constraints())

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

--------------------------------------------------------------------------------
-- Flow solver.  FE Darcy: BiCGStab + GMG (RAP, Uzawa with line Gauss-Seidel)
-- as code/flow/nephron_flow_fourth_pad100.lua.  FV1 Darcy: GMRES(100) + GMG
-- (Uzawa with ILU), 2026-11-24.
--------------------------------------------------------------------------------

if darcyScheme == "FE" then
    flowRegularInverse = LineGaussSeidel(approxSpace)
    flowRegularInverse:set_num_steps(1, 1, 1, 1, 1, 1)
else
    flowRegularInverse = ILU()
    flowRegularInverse:set_damp(1.0)
    flowRegularInverse:set_inversion_eps(1.0e-30)
end
flowPressureInverse = ILU()
flowPressureInverse:set_damp(1.0)
flowPressureInverse:set_inversion_eps(1.0e-30)

uzawa = UzawaBase({"p_stokes"})
uzawa:set_forward_iter(flowRegularInverse)
uzawa:set_schur_iter(flowPressureInverse)
uzawa:set_backward_iter(flowRegularInverse)

schurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
schurMassEq:set_mass_scale(1.0)
schurMassEq:set_diffusion(0.0)
schurMassDisc = DomainDiscretization(approxSpace)
schurMassDisc:add(schurMassEq)
schurMassDisc:add(OneSideP1Constraints())
schurMassOp = AssembledLinearOperator(schurMassDisc)
uzawa:set_schur_operator_update(schurMassOp, 1.0 / kinematicViscosity)

-- No direct (sparse LU) solver anywhere, not even on the coarsest grid: the
-- target is a very large distributed system.  The GMG coarse level is solved
-- iteratively, in parallel (not gathered), to a tight tolerance so the outer
-- (non-flexible) Krylov method sees an essentially fixed preconditioner.
-- On a run with numRefs = 0 the multigrid has a single level, so this
-- coarse solve is the whole preconditioner.
coarseSolverReduction = 1.0e-10
coarseSolverMaxIterations = 2000
function IterativeCoarseSolver()
    return {
        type = "bicgstab",
        precond = {type = "ilu"},
        convCheck = {type = "standard", iterations = coarseSolverMaxIterations,
                     absolute = 1.0e-14, reduction = coarseSolverReduction,
                     verbose = false}
    }
end

-- Known limit (2026-11-26, best mesh): at r2 this coarse solve of the
-- Galerkin coarse Stokes--Darcy operator does not reach 1e-10 in 2000
-- iterations and the steady flow solve aborts; accepting the inexact coarse
-- iterate instead makes the outer BiCGStab diverge.
flowCoarseSolverChoice = "superlu"
util.solver.defaults.approxSpace = approxSpace
if darcyScheme == "FE" then flowSolver = BiCGStab() else flowSolver = GMRES(100) end
flowSolver:set_preconditioner(util.solver.CreatePreconditioner({
    type = "gmg", approxSpace = approxSpace, discretization = domainDisc,
    rap = true, smoother = uzawa, cycle = "V", preSmooth = 2, postSmooth = 2,
    -- TEMPORARY (user, 2026-11-28): SuperLU on level 0 only, as the flow-only
    -- scripts, until an iterative coarse solve works at r2 (["superlu",
    -- IterativeCoarseSolver()]).  The transport block stays fully iterative.
    baseSolver = flowCoarseSolverChoice, baseLevel = 0,
    gatheredBaseSolverIfAmbiguous = true,
    transfer = {type = "std", restrictionDamp = 1.0, prolongationDamp = 1.0,
                enableP1LagrangeOptimization = true}
}))
flowSolver:set_convergence_check(util.solver.CreateConvCheck({
    type = "standard", iterations = 400, absolute = 1.0e-14,
    reduction = 1.0e-8, verbose = true
}))

flowMatrix = AssembledLinearOperator(domainDisc)
solution = GridFunction(approxSpace)
rightHandSide = solution:clone()
solution:set(0.0)
domainDisc:adjust_solution(solution)
solution:set_consistent_storage_type()
domainDisc:assemble_linear(flowMatrix, rightHandSide)

flowInitialDefectVector = solution:clone()
domainDisc:assemble_defect(flowInitialDefectVector, solution)
flowInitialDefect = VecNorm(flowInitialDefectVector)
flowInitialDefectVector = nil
print("Solving steady flow...")
flowSolveStart = GetClockS()
flowSolver:init(flowMatrix, solution)
if not flowSolver:apply(solution, rightHandSide) then
    print("Flow solver failed")
    exit()
end
flowSolveTime = GetClockS() - flowSolveStart
solution:set_consistent_storage_type()

-- Exact algebraic interface conservation (rows captured before the merge).
interfaceDefect = solution:clone()
pressureTrace:capture_next_interface_defect()
domainDisc:assemble_defect(interfaceDefect, solution)
interfaceSums = pressureTrace:interface_defect_sums()
stokesInterfaceFlux = interfaceSums[1]
darcyInterfaceFlux = interfaceSums[2]
interfaceMismatchPercent = 100.0 * math.abs(interfaceSums[3]) /
    math.max(math.abs(stokesInterfaceFlux), math.abs(darcyInterfaceFlux), 1e-30)

-- Darcy water flux leaving through the outer boundaries (reconstructed).
darcyMembraneVelocity = DarcyVelocityLinker()
darcyMembraneVelocity:set_permeability(membranePermeability * pressureScale)
darcyMembraneVelocity:set_viscosity(kinematicViscosity)
darcyMembraneVelocity:set_density(1.0)
darcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
darcyMembraneVelocity:set_pressure_gradient(
    ExplicitGridFunctionGradient(solution, "p_darcy_scaled"))
darcyInterVelocity = DarcyVelocityLinker()
darcyInterVelocity:set_permeability(interPermeability * pressureScale)
darcyInterVelocity:set_viscosity(kinematicViscosity)
darcyInterVelocity:set_density(1.0)
darcyInterVelocity:set_gravity(ConstUserVector(0.0))
darcyInterVelocity:set_pressure_gradient(
    ExplicitGridFunctionGradient(solution, "p_darcy_scaled"))

outerWallWater = IntegralNormalComponentOnManifold(
    darcyInterVelocity, solution, "OuterWall", "Inter")
membraneWallWater = IntegralNormalComponentOnManifold(
    darcyMembraneVelocity, solution, "MembraneWall", "Membrane")

print(string.format("Flow solver time: %.3f s", flowSolveTime))
print(string.format(
    "Assembled interface rows [mm^3/s]: Stokes=%.6e Darcy=%.6e " ..
    "(mismatch %.3g%%); reabsorbed %.4f%% of Q",
    stokesInterfaceFlux, darcyInterfaceFlux, interfaceMismatchPercent,
    100.0 * math.abs(darcyInterfaceFlux) / inletFlowRate))
print(string.format(
    "Darcy water out [mm^3/s]: OuterWall=%.6e MembraneWall=%.6e",
    outerWallWater, membraneWallWater))
if interfaceMismatchPercent > interfaceConservationTolerancePercent then
    print("Interface conservation failed")
    exit()
end
print("Steady flow stage complete")

--------------------------------------------------------------------------------
-- Transient transport discretization (block solve per step)
--------------------------------------------------------------------------------

-- Convective velocities are the face fluxes of the flow discretization
-- (Stokes velocity_ip, Darcy gradient exports), used as one-way
-- coefficients: concentrations do not feed back on the flow, so the
-- Jacobian is block lower-triangular and the flow-then-transport block solve
-- per step is exact.  velocity_ip is why flow and transport share one space.
stokesEq:set_velocity_ip_explicit(true)
-- All flow element discs are time-dependent (Stokes keeps its velocity time
-- derivative; Darcy has zero mass).  Do NOT mark only stokesEq stationary:
-- backward Euler scales the stiffness of time-dependent discs by dt, so a
-- stationary Stokes block next to time-dependent tractions, inflow and the
-- merged Darcy rows makes the flow defect inconsistent (5e3 at the steady
-- solution, 2026-11-25).  Starting from the steady flow the time derivative
-- is zero, so the flow stays quasi-steady; the same setup is the transient
-- flow model.
lumenVelocity = stokesEq:velocity_ip()

function MakeDarcyTransportVelocity(eq, permeability)
    local velocity = DarcyVelocityLinker()
    velocity:set_permeability(permeability * pressureScale)
    velocity:set_viscosity(kinematicViscosity)
    velocity:set_density(1.0)
    velocity:set_gravity(ConstUserVector(0.0))
    velocity:set_pressure_gradient(eq:gradient())
    velocity:set_derivative_mask(1)
    return velocity
end
membraneVelocity = MakeDarcyTransportVelocity(darcyMembraneEq, membranePermeability)
interVelocity = MakeDarcyTransportVelocity(darcyInterEq, interPermeability)
compartmentVelocity = {lumen = lumenVelocity, membrane = membraneVelocity,
                       inter = interVelocity}
-- NoNormalStressOutflow puts the standard nodal (P1) velocity interpolation
-- into the Outlet continuity rows; the solute outflow must use the same.
-- The coupled solution (flow and concentrations); initialised from the
-- converged steady flow, concentrations are set below.
transportSolution = solution:clone()
outletVelocity = ExplicitGridFunctionVector(transportSolution, "u,v,w")

-- Element discs that export the convective velocities.
function AddFlowExporters(disc)
    disc:add(stokesEq)
    disc:add(darcyMembraneEq)
    disc:add(darcyInterEq)
end

-- The complete Stokes--Darcy system with all its boundary conditions.
function AddFlowSystem(disc)
    AddFlowExporters(disc)
    disc:add(stokesInlet)
    disc:add(stokesOpenBoundary)
    disc:add(stokesRim)
    disc:add(stokesBJS)
    disc:add(pressureTrace)
    disc:add(darcyOuterPressure)
    disc:add(inactiveDarcy)
    disc:add(inactiveStokes)
    for _, component in ipairs({
        {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
        {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
        {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
    }) do
        local traction = NeumannBoundary(component[1], "fv1")
        traction:add(darcyPressure * component[2], "Apical", "Lumen")
        disc:add(traction)
        local outlet = NeumannBoundary(component[1], "fv1")
        outlet:add(outletKinematicPressure * component[2], "Outlet", "Lumen")
        disc:add(outlet)
    end
end

-- Convective outflow c (u.n) on open boundaries.
function MakeOutflow(field, velocityData, valueData, boundary, volume)
    local flux = ScaleAddLinkerVector()
    flux:add(valueData, velocityData)
    local neumann = NeumannBoundary(field, "fv1")
    neumann:add(flux, boundary, volume)
    return neumann
end

-- Solute influx with the inlet water.  NavierStokesInflow adds the water
-- inflow to the continuity rows as an FV1 Neumann flux of InletVelocity.n on
-- every Inlet face, including faces of the non-Dirichlet rim boxes
-- (TripleJunction).  The identical flux times c_inlet gives those boxes their
-- solute; without it they receive solute-free water and dilute.  Rows of the
-- Dirichlet Inlet vertices are overwritten by the concentration boundary.
inletVelocityData = LuaUserVector("InletVelocity")
function MakeInletInflow(field, concentration)
    local flux = ScaleAddLinkerVector()
    flux:add(ConstUserNumber(concentration), inletVelocityData)
    local neumann = NeumannBoundary(field, "fv1")
    neumann:add(flux, "Inlet", "Lumen")
    return neumann
end

-- Block solve per step.  Transport block: the flow element discs are kept
-- only to export the velocities; every flow DoF keeps its current value (the
-- flow just solved), so its rows are identities and the block contains the
-- transport rows only.
function HoldFields(boundary, fields, subsets)
    for _, field in ipairs(fields) do boundary:add(field, subsets) end
end
flowHold = DirichletBoundary()
HoldFields(flowHold, {"u", "v", "w"}, velocitySubsets)
HoldFields(flowHold, {"p_stokes"}, stokesPressureSubsets)
HoldFields(flowHold, {"p_darcy_scaled"}, allSubsets)

-- Passive leak source -> target.  Linear: J = P (c_s - c_t).  GHK at the
-- fixed potentials (v = phi_t - phi_s, k = zF/RT):
--   J = -P k v (c_s - c_t e^{kv}) / (1 - e^{kv}) = A (c_s - e^{kv} c_t),
-- which is UG4's Leak GHK formula; with fixed v it is linear in c, so a
-- Leak with permeability A and target input scaled by e^{kv} is exact.
function MakeLeak(sourceField, targetField, permeability, p,
                  sourceCompartment, targetCompartment)
    local leak = Leak({sourceField, targetField})
    local kv = p.valency * faradayConstant / (gasConstant * temperature) *
        (compartmentPotential[targetCompartment] -
         compartmentPotential[sourceCompartment])
    if not p.ghk or math.abs(kv) < 1.0e-12 then
        leak:set_permeability(permeability)
    else
        local ghkPermeability = -permeability * kv / (1.0 - math.exp(kv))
        leak:set_permeability(ghkPermeability)
        leak:set_scale_input(1, math.exp(kv))
        print(string.format("GHK leak %s %s->%s: kv=%.4f, A=%.6g mm/s " ..
            "(P=%.6g), target scale e^kv=%.6g", p.label, sourceCompartment,
            targetCompartment, kv, ghkPermeability, permeability, math.exp(kv)))
    end
    return leak
end

-- Solvent drag: J_s = (1 - sigma) J_v c_up across Apical (lumen -> cell) and
-- Basolateral (cell -> Inter).  J_v [mm/s] per node is the consistent FE
-- boundary flux of the Membrane Darcy block (its Apical sum is the assembled
-- interface water flux) divided by the node's boundary area; c_up is the
-- upstream concentration (positivity preserving).  Forward (J_v > 0) and
-- backward (J_v < 0) parts are separate transporters.  J_v is a nodal field
-- stored in the p_darcy_scaled component of auxiliary grid functions and
-- read at the vertex (MembraneTransportFV1 evaluates its density at the
-- boundary-face corner = node) by VertexGridFunctionNumberData
-- (neuro_collection, 2026-11-27): an exact vertex-position lookup.  UG4's
-- GlobalGridFunctionNumberData is unusable here: its edge search tests only
-- x (drag 20 % low) and its volume search fails at curved hex corners.
waterFluxEq = DarcyEquation("p_darcy_scaled", "Membrane")
waterFluxEq:set_mass_scale(0.0)
waterFluxEq:set_diffusion(scaledMembraneMobility)
waterFluxDisc = DomainDiscretization(approxSpace)
waterFluxDisc:add(waterFluxEq)

function NodeIndicator(field, subsets)
    local g = solution:clone()
    g:set(0.0)
    Interpolate(1.0, g, field, subsets)
    g:set_consistent_storage_type()
    return g
end

-- Boundary area per node [mm^2] (box area of the FV1 boundary faces).
function NodalBoundaryArea(boundary, indicator)
    local neumann = NeumannBoundary("p_darcy_scaled", "fv1")
    neumann:add(1.0, boundary, "Membrane")
    local disc = DomainDiscretization(approxSpace)
    disc:add(neumann)
    local area = solution:clone()
    area:set(0.0)
    disc:assemble_defect(area, solution)
    if VecProd(area, indicator) < 0.0 then VecScaleAssign(area, -1.0, area) end
    local total = VecProd(area, indicator)
    area:enforce_consistent_type()
    return area, total
end

waterIndicator = {
    -- the lumen rim (TripleJunction) nodes carry part of the Apical water
    apical = NodeIndicator("p_darcy_scaled", "Apical,TripleJunction"),
    basolateral = NodeIndicator("p_darcy_scaled", "Basolateral"),
}
waterArea, waterAreaTotal = {}, {}
for _, side in ipairs({"apical", "basolateral"}) do
    waterArea[side], waterAreaTotal[side] = NodalBoundaryArea(
        side == "apical" and "Apical" or "Basolateral", waterIndicator[side])
end
waterOnes = solution:clone()
waterOnes:set(1.0)
waterOnes:set_consistent_storage_type()
waterWork = {solution:clone(), solution:clone(), solution:clone()}
-- forward / backward flux density fields [mm/s] and their face data
waterFluxDensity = {}
waterFluxData = {}
for _, side in ipairs({"apical", "basolateral"}) do
    waterFluxDensity[side] = {forward = solution:clone(), backward = solution:clone()}
    waterFluxData[side] = {}
    for _, direction in ipairs({"forward", "backward"}) do
        waterFluxDensity[side][direction]:set(0.0)
        waterFluxData[side][direction] = VertexGridFunctionNumberData(
            waterFluxDensity[side][direction], "p_darcy_scaled")
    end
end
waterFluxResidual = solution:clone()
waterFlux = {}   -- [mm^3/s] apical: lumen -> cell, basolateral: cell -> Inter

function UpdateWaterFlux(currentSolution)
    local R = waterFluxResidual
    R:set(0.0)
    -- R_i = water entering the Membrane at node i (consistent FE flux)
    waterFluxDisc:assemble_defect(R, currentSolution)
    waterFlux.apical = VecProd(R, waterIndicator.apical)
    waterFlux.basolateral = -VecProd(R, waterIndicator.basolateral)
    R:enforce_consistent_type()
    local denom, logv, absJ = waterWork[1], waterWork[2], waterWork[3]
    for _, side in ipairs({"apical", "basolateral"}) do
        local ind = waterIndicator[side]
        local J = waterFluxDensity[side].forward
        -- 1 / area at boundary nodes (1 elsewhere)
        VecScaleAdd2(denom, 1.0, waterOnes, -1.0, ind)
        VecScaleAdd2(denom, 1.0, denom, 1.0, waterArea[side])
        VecLog(logv, denom)
        VecScaleAssign(logv, -1.0, logv)
        VecExp(denom, logv)
        -- signed J_v = (+/-) R / area on this side's nodes
        VecHadamardProd(J, R, denom)
        VecHadamardProd(J, J, ind)
        if side == "basolateral" then VecScaleAssign(J, -1.0, J) end
        -- |J_v| = sqrt(J_v^2), then forward / backward parts
        VecHadamardProd(absJ, J, J)
        VecScaleAdd2(absJ, 1.0, absJ, 1.0e-300, waterOnes)
        VecLog(logv, absJ)
        VecScaleAssign(logv, 0.5, logv)
        VecExp(absJ, logv)
        local backward = waterFluxDensity[side].backward
        VecScaleAdd2(backward, 0.5, absJ, -0.5, J)
        VecScaleAdd2(J, 0.5, absJ, 0.5, J)
        J:set_consistent_storage_type()
        backward:set_consistent_storage_type()
    end
end
UpdateWaterFlux(transportSolution)
print(string.format(
    "Solvent drag water flux [mm^3/s]: Apical %.6e (assembled interface %.6e), " ..
    "Basolateral %.6e; areas [mm^2]: Apical %.6e, Basolateral %.6e",
    waterFlux.apical, darcyInterfaceFlux, waterFlux.basolateral,
    waterAreaTotal.apical, waterAreaTotal.basolateral))

-- Drag transporter: Leak with the target input scaled to 0, so the flux is
-- (1 - sigma) c_source times the density J_v (forward or backward part).
function MakeDrag(sourceField, targetField, sigma, boundary, densityData)
    local drag = Leak({sourceField, targetField})
    drag:set_permeability(1.0 - sigma)
    drag:set_scale_input(1, 0.0)
    local disc = MembraneTransportFV1(boundary, drag)
    disc:set_density_function(densityData)
    return disc, drag
end

transportDisc = DomainDiscretization(approxSpace)
AddFlowExporters(transportDisc)
transportDisc:add(flowHold)
transportEquations = {}
for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    local eqs = {}
    for _, compartment in ipairs(compartmentOrder) do
        local eq = ConvectionDiffusionFV1(fields[compartment],
                                          compartmentVolume[compartment])
        eq:set_diffusion(p.diffusion)
        eq:set_velocity(compartmentVelocity[compartment])
        eq:set_upwind(FullUpwind())
        transportDisc:add(eq)
        eqs[compartment] = eq
    end

    eqs.membraneDiscs = {}
    if p.apicalPermeability > 0.0 then
        local apicalLeak = MakeLeak(fields.lumen, fields.membrane,
                                    p.apicalPermeability, p,
                                    "lumen", "membrane")
        local apical = MembraneTransportFV1("Apical", apicalLeak)
        apical:set_density_function(1.0)
        transportDisc:add(apical)
        eqs.membraneDiscs[#eqs.membraneDiscs + 1] = apical
    end
    if p.basolateralPermeability > 0.0 then
        local basolateralLeak = MakeLeak(fields.membrane, fields.inter,
                                         p.basolateralPermeability, p,
                                         "membrane", "inter")
        local basolateral = MembraneTransportFV1("Basolateral", basolateralLeak)
        basolateral:set_density_function(1.0)
        transportDisc:add(basolateral)
        eqs.membraneDiscs[#eqs.membraneDiscs + 1] = basolateral
    end

    -- solvent drag (forward: with the water, backward: water flowing back)
    eqs.dragDiscs = {apical = {}, basolateral = {}}
    eqs.dragObjects = {}
    for _, d in ipairs({
        {"apical", "Apical", fields.lumen, fields.membrane, p.reflectionApical},
        {"basolateral", "Basolateral", fields.membrane, fields.inter,
         p.reflectionBasolateral}}) do
        local side, boundary, inner, outer, sigma = d[1], d[2], d[3], d[4], d[5]
        if sigma < 1.0 then
            local fwd, fwdT = MakeDrag(inner, outer, sigma, boundary,
                                       waterFluxData[side].forward)
            local bwd, bwdT = MakeDrag(outer, inner, sigma, boundary,
                                       waterFluxData[side].backward)
            for _, disc in ipairs({fwd, bwd}) do
                transportDisc:add(disc)
                eqs.membraneDiscs[#eqs.membraneDiscs + 1] = disc
            end
            eqs.dragDiscs[side] = {fwd, bwd}
            eqs.dragObjects[#eqs.dragObjects + 1] = fwdT
            eqs.dragObjects[#eqs.dragObjects + 1] = bwdT
        end
    end

    eqs.neumannDiscs = {
        MakeInletInflow(fields.lumen, p.inlet),
        MakeOutflow(fields.lumen, outletVelocity, eqs.lumen:value(),
                    "Outlet", "Lumen")}
    for _, neumann in ipairs(eqs.neumannDiscs) do transportDisc:add(neumann) end
    transportEquations[soluteKey] = eqs
end
pumpDiscs = {}
if pumpEnabled then
    naKPump = NaKPump({transportFields.na.membrane, transportFields.na.inter,
                       transportFields.k.membrane, transportFields.k.inter})
    naKPump:set_max_flux(pumpMaxCycleFlux)
    naKPump:set_k_na(pumpKNa)
    naKPump:set_k_k(pumpKK)
    local pump = MembraneTransportFV1("Basolateral", naKPump)
    pump:set_density_function(1.0)
    transportDisc:add(pump)
    pumpDiscs[1] = pump
end

-- Inlet: inlet concentration.  OuterWall: fixed bath (outerWall).  The
-- Darcy pressure is Dirichlet on OuterWall, so the water leaving
-- there is the internal FV1 face flux into those boxes; a fixed bath
-- concentration is the consistent (and natural) transport condition.
concentrationBoundary = DirichletBoundary()
for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    concentrationBoundary:add(p.inlet, transportFields[soluteKey].lumen, "Inlet")
    concentrationBoundary:add(p.outerWall, transportFields[soluteKey].inter,
                              "OuterWall")
end
transportDisc:add(concentrationBoundary)
transportDisc:add(OneSideP1Constraints())

-- Flow block: the full Stokes--Darcy system (stationary, quasi-steady flow)
-- with every concentration held.  Its Jacobian equals the steady-flow
-- matrix, so the validated flow GMG (flowSolver) is reused.
concentrationHold = DirichletBoundary()
for _, soluteKey in ipairs(soluteOrder) do
    for _, compartment in ipairs(compartmentOrder) do
        concentrationHold:add(transportFields[soluteKey][compartment],
                              transportFieldSubsets[compartment])
    end
end
flowStepDisc = DomainDiscretization(approxSpace)
AddFlowSystem(flowStepDisc)
flowStepDisc:add(concentrationHold)
flowStepDisc:add(OneSideP1Constraints())

for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    Interpolate(p.inlet, transportSolution, fields.lumen,
                transportFieldSubsets.lumen)
    Interpolate(p.cellInitial, transportSolution, fields.membrane,
                transportFieldSubsets.membrane)
    Interpolate(p.interInitial, transportSolution, fields.inter,
                transportFieldSubsets.inter)
end
transportDisc:adjust_solution(transportSolution)
transportSolution:set_consistent_storage_type()

--------------------------------------------------------------------------------
-- Discrete mass-balance diagnostics
--------------------------------------------------------------------------------
-- Lumped FV1 mass M = sum_i |box_i| c_i (exactly the transient mass term).
-- R = spatial residual of all transport element discs (no concentration
-- Dirichlet), evaluated exactly as in the time step: R_i = D_i/dt at the
-- Dirichlet nodes, with D the backward-Euler defect of these discs for the
-- step just solved (the flow is time-dependent, so its FV1 stabilization and
-- hence velocity_ip depend on dt and the old solution; a separately assembled
-- stationary residual would use a different velocity).
-- Backward Euler gives |box_i|(c_i-c_i^old)/dt + R_i = 0 at every
-- non-Dirichlet node, and c is fixed at Dirichlet nodes, hence
--   dM/dt = -sum_all R + sum_Dirichlet R.
-- Internal, diffusive and two-sided membrane fluxes cancel in sum_all R,
-- which leaves the Neumann rows.  Hence, exactly,
--   inflow  = -sum(inlet Neumann rows) + sum_{Inlet nodes} R
--   outlet  =  sum(outlet Neumann rows)
--   outer   = -sum_{OuterWall nodes} R      (net solute leaving to the bath)
--   dM/dt   = inflow - outlet - outer.
-- Any one-sided (non-cancelling) term appears as a balance error.

massDisc = DomainDiscretization(approxSpace)
stationaryDisc = DomainDiscretization(approxSpace)
inletNeumannDisc = DomainDiscretization(approxSpace)
outletNeumannDisc = DomainDiscretization(approxSpace)
AddFlowExporters(stationaryDisc)
for _, soluteKey in ipairs(soluteOrder) do
    local p = soluteParameters[soluteKey]
    local fields = transportFields[soluteKey]
    for _, compartment in ipairs(compartmentOrder) do
        local m = ConvectionDiffusionFV1(fields[compartment],
                                         compartmentVolume[compartment])
        m:set_diffusion(0.0)
        m:set_reaction_rate(1.0)
        massDisc:add(m)
    end
    inletNeumannDisc:add(MakeInletInflow(fields.lumen, p.inlet))
    local outletValue = ConvectionDiffusionFV1(fields.lumen, "Lumen")
    outletValue:set_diffusion(0.0)  -- zero rows; only provides value()
    outletNeumannDisc:add(outletValue)
    outletNeumannDisc:add(MakeOutflow(fields.lumen, outletVelocity,
                                      outletValue:value(), "Outlet", "Lumen"))
end
-- the stationary residual uses the transport element discs themselves
for _, pump in ipairs(pumpDiscs) do stationaryDisc:add(pump) end
for _, soluteKey in ipairs(soluteOrder) do
    local eqs = transportEquations[soluteKey]
    for _, compartment in ipairs(compartmentOrder) do
        stationaryDisc:add(eqs[compartment])
    end
    for _, membrane in ipairs(eqs.membraneDiscs) do
        stationaryDisc:add(membrane)
    end
    for _, neumann in ipairs(eqs.neumannDiscs) do
        stationaryDisc:add(neumann)
    end
end

function Indicator(field, subsets)
    local g = solution:clone()
    g:set(0.0)
    Interpolate(1.0, g, field, subsets)
    g:set_consistent_storage_type()
    return g
end
indicators = {}
for _, soluteKey in ipairs(soluteOrder) do
    local fields = transportFields[soluteKey]
    local all = solution:clone()
    all:set(0.0)
    for _, compartment in ipairs(compartmentOrder) do
        Interpolate(1.0, all, fields[compartment],
                    transportFieldSubsets[compartment])
    end
    all:set_consistent_storage_type()
    indicators[soluteKey] = {
        all = all,
        inlet = Indicator(fields.lumen, "Inlet"),
        outer = Indicator(fields.inter, "OuterWall"),
        lumen = Indicator(fields.lumen, transportFieldSubsets.lumen),
        membrane = Indicator(fields.membrane, transportFieldSubsets.membrane),
        inter = Indicator(fields.inter, transportFieldSubsets.inter),
    }
end

-- Solvent-drag-only discretizations, for the drag part of the transfers.
dragBudgetDisc = {}
for _, soluteKey in ipairs(soluteOrder) do
    dragBudgetDisc[soluteKey] = {}
    for _, side in ipairs({"apical", "basolateral"}) do
        local discs = transportEquations[soluteKey].dragDiscs[side]
        if #discs > 0 then
            local d = DomainDiscretization(approxSpace)
            for _, disc in ipairs(discs) do d:add(disc) end
            dragBudgetDisc[soluteKey][side] = d
        end
    end
end

diagDefect = solution:clone()

budgetTimeDisc = ThetaTimeStep(stationaryDisc, 1.0)

-- oldSeries/stepDt: the time series and step of the step just solved (before
-- it is pushed); nil at t = 0 (mass only).
function SoluteBudget(currentSolution, oldSeries, stepDt)
    local budget = {}
    diagDefect:set(0.0)
    massDisc:assemble_defect(diagDefect, currentSolution)
    for _, soluteKey in ipairs(soluteOrder) do
        local ind = indicators[soluteKey]
        budget[soluteKey] = {mass = VecProd(diagDefect, ind.all),
                             massLumen = VecProd(diagDefect, ind.lumen),
                             massMembrane = VecProd(diagDefect, ind.membrane),
                             massInter = VecProd(diagDefect, ind.inter)}
    end
    if oldSeries == nil then return budget end
    -- solvent-drag rates (net, source side of each membrane)
    for _, soluteKey in ipairs(soluteOrder) do
        local b, ind = budget[soluteKey], indicators[soluteKey]
        b.dragApical, b.dragBasolateral = 0.0, 0.0
        for side, d in pairs(dragBudgetDisc[soluteKey]) do
            diagDefect:set(0.0)
            d:assemble_defect(diagDefect, currentSolution)
            if side == "apical" then
                b.dragApical = VecProd(diagDefect, ind.lumen)
            else
                b.dragBasolateral = VecProd(diagDefect, ind.membrane)
            end
        end
    end
    diagDefect:set(0.0)
    budgetTimeDisc:prepare_step(oldSeries, stepDt)
    budgetTimeDisc:assemble_defect(diagDefect, currentSolution)
    for _, soluteKey in ipairs(soluteOrder) do
        local b = budget[soluteKey]
        b.inletDirichlet = VecProd(diagDefect, indicators[soluteKey].inlet) / stepDt
        b.outer = -VecProd(diagDefect, indicators[soluteKey].outer) / stepDt
    end
    diagDefect:set(0.0)
    inletNeumannDisc:assemble_defect(diagDefect, currentSolution)
    for _, soluteKey in ipairs(soluteOrder) do
        local b = budget[soluteKey]
        b.inflow = b.inletDirichlet - VecProd(diagDefect, indicators[soluteKey].lumen)
    end
    diagDefect:set(0.0)
    outletNeumannDisc:assemble_defect(diagDefect, currentSolution)
    for _, soluteKey in ipairs(soluteOrder) do
        budget[soluteKey].outlet = VecProd(diagDefect, indicators[soluteKey].lumen)
    end
    return budget
end

-- Positivity at vertex DoFs.  (A convex-hull bound does not hold: with no
-- solvent drag, water crossing Apical/Basolateral leaves its solute behind,
-- which concentrates the source side; ranges are checked from the VTU output.)
function CheckBounds(currentSolution)
    local ok = true
    for _, soluteKey in ipairs(soluteOrder) do
        for _, compartment in ipairs(compartmentOrder) do
            local field = transportFields[soluteKey][compartment]
            if not CheckGFValuesWithinBounds(currentSolution, field, 0.0, 1.0e6) then
                print("POSITIVITY FAILED: " .. field)
                ok = false
            end
        end
    end
    return ok
end

previousBudget = SoluteBudget(transportSolution)
maxBalanceErrorPercent = 0.0
for _, soluteKey in ipairs(soluteOrder) do
    local b = previousBudget[soluteKey]
    print(string.format("BUDGET t=0 %s: mass=%.10e (mM mm^3)",
                        soluteParameters[soluteKey].label, b.mass))
end

--------------------------------------------------------------------------------
-- Block solvers
--------------------------------------------------------------------------------
-- A monolithic GMG on all 17 fields stalled at r1 (coarse Uzawa solve did
-- not converge, flow rows erratic; 2026-11-25/26).  Each step therefore
-- solves (1) the flow block with the validated flow GMG and (2) the
-- transport block with its own GMG.  Both blocks are linear: the transport
-- Jacobian is assembled once (rebuilt only if the flow changes) and each
-- step performs defect correction u <- u - J^{-1} d(u) with the exact defect.

timeDisc = ThetaTimeStep(transportDisc, 1.0)

-- Transport block: BiCGStab + GMG (Galerkin coarse operators, ILU smoother,
-- iterative BiCGStab/ILU coarse solve; no direct solver).  On the 10 um
-- shell mesh this converges in 6 (r1) and 8-10 (r2) iterations for the
-- diffusion+leak test (2026-11-26).
util.solver.defaults.approxSpace = approxSpace
transportLinearSolver = util.solver.CreateSolver({
    type = "bicgstab",
    precond = {
        type = "gmg", approxSpace = approxSpace, discretization = timeDisc,
        baseLevel = 0, rap = true, cycle = "V", preSmooth = 3, postSmooth = 3,
        smoother = {type = "ilu", overlap = true},
        transfer = "std",
        baseSolver = IterativeCoarseSolver(),
    },
    convCheck = {type = "standard", iterations = 200, absolute = 1.0e-14,
                 reduction = 1.0e-10, verbose = false},
})
transportDefectAbsolute = 1.0e-10
-- The Na/K pump makes the transport block nonlinear.  Each step iterates
-- u <- u - J^{-1} d(u) with the exact defect and a reused Jacobian
-- (chord/Newton); J and the GMG are rebuilt at the current state when the
-- previous step needed more than jacobianRebuildCorrections corrections.
transportMaxCorrections = 10
jacobianRebuildCorrections = 2

-- Flow block per step: solved when its defect exceeds flowStepTolerance
-- (relative 1e-8 of the steady-flow initial defect, as the steady solve).
flowStepConv = ConvCheck(400, 1.0e-14, 1.0e-8, false)
flowStepTolerance = 1.0e-8 * flowInitialDefect

--------------------------------------------------------------------------------
-- Compartment output
--------------------------------------------------------------------------------

outputTimes = {}
lastOutputStep = 0
compartmentOutputs = {}
for _, compartment in ipairs(compartmentOrder) do
    local writer = VTKOutput()
    writer:clear_selection()
    for _, soluteKey in ipairs(soluteOrder) do
        writer:select(transportFields[soluteKey][compartment],
                      soluteParameters[soluteKey].label)
    end
    writer:select_element(compartmentVelocity[compartment], "velocity_mm_s")
    compartmentOutputs[#compartmentOutputs + 1] = {
        name = compartment, subset = compartmentVolume[compartment],
        writer = writer}
end

function WriteSerialPVTU(prefix, step, time)
    if ProcRank() ~= 0 or NumProcs() ~= 1 then return end
    local file = assert(io.open(string.format("%s_t%04d.pvtu", prefix, step), "w"))
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n<VTKFile type="PUnstructuredGrid" ',
               'version="0.1" byte_order="LittleEndian">\n')
    file:write(string.format('  <Time timestep="%.17g"/>\n', time))
    file:write('  <PUnstructuredGrid GhostLevel="0">\n    <PPoints>',
               '<PDataArray type="Float32" NumberOfComponents="3"/></PPoints>\n',
               '    <PPointData>\n')
    for _, soluteKey in ipairs(soluteOrder) do
        file:write(string.format('      <PDataArray type="Float32" Name="%s" ' ..
            'NumberOfComponents="1"/>\n', soluteParameters[soluteKey].label))
    end
    file:write('    </PPointData>\n    <PCellData><PDataArray type="Float32" ',
               'Name="velocity_mm_s" NumberOfComponents="3"/></PCellData>\n')
    file:write(string.format('    <Piece Source="%s_t%04d.vtu"/>\n', basename, step))
    file:write('  </PUnstructuredGrid>\n</VTKFile>\n')
    file:close()
end

function WriteCompartmentPVD(prefix)
    if ProcRank() ~= 0 then return end
    local file = assert(io.open(prefix .. ".pvd", "w"))
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" ',
               'byte_order="LittleEndian">\n  <Collection>\n')
    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            file:write(string.format('    <DataSet timestep="%.17g" group="" ' ..
                'part="0" file="%s_t%04d.pvtu"/>\n', outputTimes[step], basename, step))
        end
    end
    file:write('  </Collection>\n</VTKFile>\n')
    file:close()
end

function WriteTransportOutput(currentSolution, step, time)
    for _, c in ipairs(compartmentOutputs) do
        local prefix = outputName .. "_" .. c.name
        c.writer:print_subsets(prefix, currentSolution, c.subset, step, time)
        WriteSerialPVTU(prefix, step, time)
    end
    outputTimes[step] = time
    lastOutputStep = math.max(lastOutputStep, step)
    for _, c in ipairs(compartmentOutputs) do
        WriteCompartmentPVD(outputName .. "_" .. c.name)
    end
end

--------------------------------------------------------------------------------
-- Time loop
--------------------------------------------------------------------------------

WriteTransportOutput(transportSolution, 0, 0.0)
CheckBounds(transportSolution)

function TransportStepCallback(step, time, currentDt, current, oldSeries)
    local budget = SoluteBudget(current, oldSeries, currentDt)
    for _, soluteKey in ipairs(soluteOrder) do
        local b = budget[soluteKey]
        local rate = (b.mass - previousBudget[soluteKey].mass) / currentDt
        local net = b.inflow - b.outlet - b.outer
        local scale = math.max(math.abs(b.inflow), math.abs(b.outlet),
                               math.abs(b.outer), 1.0e-30)
        local errorPercent = 100.0 * math.abs(rate - net) / scale
        maxBalanceErrorPercent = math.max(maxBalanceErrorPercent, errorPercent)
        print(string.format(
            "BUDGET step=%d t=%.6g %s: mass=%.10e dM/dt=%.6e in=%.6e " ..
            "outlet=%.6e outer=%.6e balance_err=%.3e%%",
            step, time, soluteParameters[soluteKey].label, b.mass, rate,
            b.inflow, b.outlet, b.outer, errorPercent))
    end
    -- reabsorption bookkeeping [pmol]: 1 mM mm^3 = 1 nmol = 1000 pmol
    local w = 1000.0 * currentDt
    cumulativeWater.apical = cumulativeWater.apical + waterFlux.apical * currentDt
    cumulativeWater.basolateral = cumulativeWater.basolateral +
        waterFlux.basolateral * currentDt
    for _, soluteKey in ipairs(soluteOrder) do
        local b, o, c = budget[soluteKey], previousBudget[soluteKey],
            cumulative[soluteKey]
        local apical = b.inflow - b.outlet -
            (b.massLumen - o.massLumen) / currentDt
        local basolateral = apical - (b.massMembrane - o.massMembrane) / currentDt
        c.inflow = c.inflow + w * b.inflow
        c.outlet = c.outlet + w * b.outlet
        c.bath = c.bath + w * b.outer
        c.apical = c.apical + w * apical
        c.apicalDrag = c.apicalDrag + w * b.dragApical
        c.basolateral = c.basolateral + w * basolateral
        c.basolateralDrag = c.basolateralDrag + w * b.dragBasolateral
        if ProcRank() == 0 then
            reabsorptionFile:write(string.format(
                "%.6g,%s,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e," ..
                "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                time, soluteParameters[soluteKey].label,
                1000.0 * b.inflow, 1000.0 * b.outlet, 1000.0 * apical,
                1000.0 * b.dragApical, 1000.0 * basolateral,
                1000.0 * b.dragBasolateral, 1000.0 * b.outer,
                c.inflow, c.outlet, c.apical, c.apicalDrag, c.basolateral,
                c.basolateralDrag, c.bath))
            reabsorptionFile:flush()
        end
    end
    if step % outputEverySteps == 0 or step == numSteps then
        PrintReabsorption(time, budget)
    end
    previousBudget = budget
    if not CheckBounds(current) then positivityFailed = true end
    if step % outputEverySteps == 0 or step == numSteps then
        WriteTransportOutput(current, step, time)
    end
end

-- Cumulative reabsorption [pmol] and water [mm^3] since t = 0.
cumulativeWater = {apical = 0.0, basolateral = 0.0}
cumulative = {}
for _, soluteKey in ipairs(soluteOrder) do
    cumulative[soluteKey] = {inflow = 0.0, outlet = 0.0, bath = 0.0, apical = 0.0,
                             apicalDrag = 0.0, basolateral = 0.0,
                             basolateralDrag = 0.0}
end
reabsorptionFileName = outputName .. "_reabsorption.csv"
if ProcRank() == 0 then
    reabsorptionFile = io.open(reabsorptionFileName, "w")
    if reabsorptionFile == nil then
        print("Cannot create " .. reabsorptionFileName)
        exit()
    end
    reabsorptionFile:write(
        "time_s,solute,rate_in_pmol_s,rate_outlet_pmol_s,rate_apical_pmol_s," ..
        "rate_apical_drag_pmol_s,rate_basolateral_pmol_s," ..
        "rate_basolateral_drag_pmol_s,rate_bath_pmol_s," ..
        "cum_in_pmol,cum_outlet_pmol,cum_apical_pmol,cum_apical_drag_pmol," ..
        "cum_basolateral_pmol,cum_basolateral_drag_pmol,cum_bath_pmol\n")
end

-- apical = net lumen -> cell (leak + drag + transport), basolateral = net
-- cell -> Inter; reabsorbed fraction = apical / inflow (cumulative).
function PrintReabsorption(time, budget)
    print(string.format(
        "REABSORPTION t=%.6g s  water: Apical %.6e mm^3 (%.4f%% of inflow), " ..
        "Basolateral %.6e mm^3", time, cumulativeWater.apical,
        100.0 * waterFlux.apical / inletFlowRate, cumulativeWater.basolateral))
    for _, soluteKey in ipairs(soluteOrder) do
        local c, b = cumulative[soluteKey], budget[soluteKey]
        print(string.format(
            "REABSORPTION t=%.6g s %-4s [pmol] in=%.6e outlet=%.6e " ..
            "apical=%.6e (drag %.6e) basolateral=%.6e (drag %.6e) bath=%.6e " ..
            "reabsorbed=%.3f%% | stored lumen=%.6e cell=%.6e inter=%.6e",
            time, soluteParameters[soluteKey].label, c.inflow, c.outlet,
            c.apical, c.apicalDrag, c.basolateral, c.basolateralDrag, c.bath,
            100.0 * c.apical / math.max(c.inflow, 1.0e-30),
            1000.0 * b.massLumen, 1000.0 * b.massMembrane, 1000.0 * b.massInter))
    end
end

solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(transportSolution:clone(), 0.0)
gridLevel = transportSolution:grid_level()
transportJacobian = AssembledLinearOperator(timeDisc, gridLevel)
transportDefect = transportSolution:clone()
transportCorrection = transportSolution:clone()

flowDefect = transportSolution:clone()
flowCorrection = transportSolution:clone()
flowSolver:set_convergence_check(flowStepConv)

-- The transport Jacobian depends on the (held) velocities; it is assembled
-- once and rebuilt whenever the flow block changes the flow.
function SetupTransportSolver()
    timeDisc:assemble_jacobian(transportJacobian, transportSolution, gridLevel)
    assert(transportLinearSolver:init(transportJacobian, transportSolution),
           "Transport linear solver init failed")
end

print("Assembling transport Jacobian and initializing GMG...")
transportStart = GetClockS()
timeDisc:prepare_step(solTimeSeries, dt)
SetupTransportSolver()
print(string.format("Transport setup time: %.3f s", GetClockS() - transportStart))

positivityFailed = false
solverFailed = false
numSteps = math.floor(endTime / dt + 0.5)
time = 0.0
print("Starting transient transport solve...")
for step = 1, numSteps do
    timeDisc:prepare_step(solTimeSeries, dt)
    timeDisc:prepare_step_elem(solTimeSeries, dt)
    local converged = true
    local log = {}

    -- (1) Flow block (stationary Stokes--Darcy, concentrations held).
    local flowStart = GetClockS()
    flowStepDisc:assemble_defect(flowDefect, transportSolution)
    local flowNorm = VecNorm(flowDefect)
    local flowLog = string.format("flow %.3e", flowNorm)
    if not (flowNorm == flowNorm) then converged = false end
    if converged and flowNorm > flowStepTolerance then
        flowCorrection:set(0.0)
        if not flowSolver:apply_return_defect(flowCorrection, flowDefect) then
            converged = false
            flowLog = flowLog .. " flow-solver-failed"
        end
        VecScaleAdd2(transportSolution, 1.0, transportSolution,
                     -1.0, flowCorrection)
        flowStepDisc:assemble_defect(flowDefect, transportSolution)
        flowLog = flowLog .. string.format(" (lin %d) -> %.3e",
                                           flowStepConv:step(), VecNorm(flowDefect))
        UpdateWaterFlux(transportSolution)
        SetupTransportSolver()
    else
        flowLog = flowLog .. " (held)"
    end
    local flowTime = GetClockS() - flowStart

    -- (2) Transport block (flow DoFs held at the new flow).
    local transportStepStart = GetClockS()
    if (previousCorrections or 0) > jacobianRebuildCorrections then
        SetupTransportSolver()
        log[#log + 1] = "[J rebuilt]"
    end
    local corrections = 0
    for correction = 0, transportMaxCorrections do
        timeDisc:assemble_defect(transportDefect, transportSolution, gridLevel)
        local norm = VecNorm(transportDefect)
        log[#log + 1] = string.format("%.3e", norm)
        if not (norm == norm) or norm > 1.0e30 then converged = false; break end
        if norm < transportDefectAbsolute then break end
        if correction == transportMaxCorrections then converged = false; break end
        transportCorrection:set(0.0)
        if not transportLinearSolver:apply_return_defect(
                transportCorrection, transportDefect) then
            converged = false
            log[#log + 1] = "linear-solver-failed"
        end
        log[#log + 1] = string.format("(lin %d)", transportLinearSolver:step())
        VecScaleAdd2(transportSolution, 1.0, transportSolution,
                     -1.0, transportCorrection)
        corrections = corrections + 1
    end
    previousCorrections = corrections
    time = solTimeSeries:time(0) + dt  -- same rounding as prepare_step
    print(string.format("STEP %d t=%.6g | %s, %.2f s | transport defect %s, %.2f s | %s",
                        step, time, flowLog, flowTime, table.concat(log, " "),
                        GetClockS() - transportStepStart,
                        converged and "converged" or "NOT CONVERGED"))
    if not converged then solverFailed = true; break end

    -- budget uses the old solution series, so evaluate it before the push
    TransportStepCallback(step, time, dt, transportSolution, solTimeSeries)

    local oldest = solTimeSeries:oldest()
    VecScaleAssign(oldest, 1.0, transportSolution)
    solTimeSeries:push_discard_oldest(oldest, time)
    timeDisc:finish_step_elem(solTimeSeries, gridLevel)
end

if ProcRank() == 0 then reabsorptionFile:close() end
print("Reabsorption time series: " .. reabsorptionFileName)
print(string.format("Transport time: %.3f s; total %.3f s",
                    GetClockS() - transportStart, GetClockS() - totalStartTime))
print(string.format("Maximum per-step mass-balance error: %.3e%%",
                    maxBalanceErrorPercent))
print("Positivity: " .. (positivityFailed and "FAILED" or "ok"))
print("Transport output prefix: " .. outputName)
if solverFailed then print("Transport solver: FAILED") end
if solverFailed or positivityFailed or maxBalanceErrorPercent > massBalanceTolerancePercent then
    print("STANDALONE_FLOW_TRANSPORT: COMPLETE WITH FAILED CHECKS")
else
    print("STANDALONE_FLOW_TRANSPORT: COMPLETE")
end
