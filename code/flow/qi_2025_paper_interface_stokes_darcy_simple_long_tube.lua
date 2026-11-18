--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Experimental reproduction of the interface formulation in Qi et al.
-- (International Journal for Numerical Methods in Biomedical Engineering,
-- 2025, e70108). This is an unvalidated straight-tube debugging control using
-- the paper's three interface conditions:
-- normal mass conservation, normal-stress balance, and tangential BJS slip.
-- It deliberately omits the artificial membrane-normal Robin impedance.  A
-- small algebraic defect or lumen balance does not validate the Darcy velocity;
-- the reconstructed porous boundary flux must also converge under refinement.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolves paths, loads UG4 utilities, and requires the discretization,
-- discretization and coarse-grid solver plugins used below.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({
    "NavierStokes", "ConvectionDiffusion", "neuro_collection", "SuperLU6"
})

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------
-- Only -numRefs and -output are terminal inputs. The mesh, physical model,
-- pressure scaling, boundary data, and GMG controls are configured here.

dim = 3

-- Physics, discretization, and solver choices are edited here. The only
-- command-line inputs are -numRefs and -output.
numRefs = util.GetParamNumber("-numRefs", 0, "number of uniform refinements")
bstokes = true                      -- true: Stokes; false: Navier--Stokes
vmax = 3.0                          -- maximum inlet velocity [mm/s]

-- Geometry and physical parameters: edit values here, not on the command line.
gridS = "straight_nephron_no_bend_plain"
gridName = "runs/" .. gridS .. ".ugx"
numPreRefs = 0
-- Mesh-validated straight lumen: z in [-0.42, 0.42] mm, centered at the
-- origin, with the Apical wall at r = 0.015 mm.
inletRadius = 1.5e-2               -- [mm]
inletCenterX = 0.0                 -- [mm]
inletCenterY = 0.0                 -- [mm]
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
membranePermeabilityM2 = 1.0e-15  -- membrane Darcy permeability [m^2]
interPermeabilityM2 = 1.0e-14     -- interstitial Darcy permeability [m^2]
bjsInterfacePermeabilityM2 = 1.0e-12 -- effective BJS slip permeability [m^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- [kg/m^3]
outerPressurePa = 0.0             -- [Pa], no imposed interstitial back-pressure
-- Symmetric algebraic scaling for the Darcy/interface block.  Scaling the
-- unknown alone makes the two monolithic coupling blocks badly unbalanced.
-- With S=R=1/sqrt(K_mem/nu), the membrane diffusion R*K*S/nu is order one and
-- the pressure-to-traction and velocity-to-flux couplings have equal scale.
darcyPressureScale = 1.0 / math.sqrt(membranePermeabilityM2 * 1.0e6 / viscosity)
darcyEquationScale = darcyPressureScale
solverSetup = "gmg"
gmgUseRAP = false                 -- rediscretize on each level; RAP is unstable across r2 vertical copies
interfaceFluxSign = 1.0
interfaceTractionSign = 1.0
maxKrylovCycles = 8                -- allow verified true-defect restarts
-- r1 has separately been verified with a level-0 coarse solve.  For the r2
-- mesh-convergence diagnostic, use that fixed r1 system as the base to isolate
-- r2 smoothing from the poor r0-to-r1 coarse correction.  This moderate direct
-- base is diagnostic only and must return to level 0 for production geometry.
gmgBaseLevel = math.min(1, numRefs)

-- Keep unverified debugging output out of Results_final.  MPI runs create a
-- .pvtu master file and one .vtu file per rank.
vtk_file_name = util.GetParam(
    "-output",
    "tmp/coupled_tube_r" .. numRefs,
    "VTK output prefix"
)

assert(numRefs >= 0, "numRefs must be non-negative")
assert(vmax >= 0.0, "vmax must be non-negative")
assert(membranePermeabilityM2 > 0.0,
    "membranePermeabilityM2 must be positive")
assert(interPermeabilityM2 > 0.0,
    "interPermeabilityM2 must be positive")
assert(bjsInterfacePermeabilityM2 > 0.0,
    "bjsInterfacePermeabilityM2 must be positive")
assert(bjsAlpha > 0.0, "bjsAlpha must be positive")
assert(density > 0.0, "density must be positive")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print(string.format(
    "Model: %s, vmax=%.6e mm/s, refinements=%d",
    bstokes and "Stokes" or "Navier--Stokes", vmax, numRefs
))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Loads the deliberately plain straight-tube control mesh, validates its
-- subsets, refines it geometrically, and creates one coupled space.  This UGX
-- has no ProjectionHandler by design; it is used to isolate the flow coupling
-- from curved-geometry projection effects.

print("Plain straight-tube control grid: " .. gridName)

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Membrane", "Inter", "Apical", "Basolateral",
     "OuterWall", "MembraneWall", "Inlet", "Outlet", "TripleJunction"}),
    "The domain is missing a required Stokes--Darcy subset")

-- Level-0 distribution triggers an MPI message-truncation bug in this build's
-- mixed approximation-space construction.  Distribute at level 1 instead:
-- r2 is then refined from an already partitioned r1 hierarchy without the
-- nearly full-grid finest-level ghosts caused by delaying until level 2.
balancer.firstDistLvl = math.min(1, numRefs)
-- Keep all load-balancing choices in this script; do not expose UG4's
-- additional balancing command-line parameters.
balancer.parametersParsed = true
print("First distribution level: " .. balancer.firstDistLvl)
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- Approximation-space fields
--------------------------------------------------------------------------------

approxSpace = ApproximationSpace(dom)
local allFlowSubsets =
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet,TripleJunction"
local stokesVelocitySubsetNames = {
    "Lumen", "Membrane", "Apical", "Basolateral", "MembraneWall",
    "Inlet", "Outlet", "TripleJunction"
}
approxSpace:add_fct(
    {"u", "v", "w"}, "Lagrange", 1, stokesVelocitySubsetNames
)
approxSpace:add_fct(
    "p_stokes", "Lagrange", 1, allFlowSubsets
)
approxSpace:add_fct(
    "p_darcy_scaled", "Lagrange", 1, allFlowSubsets
)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Monolithic Stokes--Darcy approximation space:")
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Defines mesh-derived inlet and interface directions first, then assembles
-- the monolithic Stokes--Darcy, normal-continuity, and BJS terms.

local membranePermeability = membranePermeabilityM2 * 1.0e6 -- [mm^2]
local interPermeability = interPermeabilityM2 * 1.0e6       -- [mm^2]
local bjsInterfacePermeability = bjsInterfacePermeabilityM2 * 1.0e6
local membraneDarcyMobility = membranePermeability / viscosity
local interDarcyMobility = interPermeability / viscosity
local scaledMembraneDarcyMobility =
    darcyEquationScale * membraneDarcyMobility * darcyPressureScale
local scaledInterDarcyMobility =
    darcyEquationScale * interDarcyMobility * darcyPressureScale
local bjsFriction =
    bjsAlpha * viscosity / math.sqrt(bjsInterfacePermeability)
local outerPressurePhysical = outerPressurePa / (density * 1.0e-6)
local outerPressureScaled = outerPressurePhysical / darcyPressureScale

print(string.format(
    "K_mem=%.6e mm^2, K_inter=%.6e mm^2, pressure scale=%.6e, Darcy row scale=%.6e",
    membranePermeability, interPermeability, darcyPressureScale,
    darcyEquationScale
))
print(string.format(
    "Scaled Darcy diffusion: membrane=%.6e, inter=%.6e; BJS friction=%.6e; paper normal Robin=0; OuterWall pressure=%.6e Pa",
    scaledMembraneDarcyMobility, scaledInterDarcyMobility,
    bjsFriction, outerPressurePa
))

function inletVelocity(x, y, z, t)
    local dx, dy = x - inletCenterX, y - inletCenterY
    local radiusSquared = dx * dx + dy * dy
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return 0.0, 0.0, vmax * profile
end

-- Stokes/Navier--Stokes equation in Lumen
--------------------------------------------------------------------------------

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
stokesEq:set_stokes(bstokes)
stokesEq:set_laplace(false) -- symmetric-gradient stress required for BJS
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_upwind("lps")
stokesEq:set_peclet_blend(false)
stokesEq:set_stabilization("fields", "cor")
stokesEq:set_pac_upwind(false)

stokesInlet = NavierStokesInflow(stokesEq)
stokesInlet:add("inletVelocity", "Inlet")
stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
-- This boundary discretization is also required on Apical: besides its
-- traction terms, it closes the Stokes continuity control volumes with the
-- actual boundary-normal velocity.  WSBC alone supplies only the tangential
-- BJS term and does not assemble that mass flux.
stokesOpenBoundary:add("Outlet,Apical")
stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

-- Qi et al., eq. (5): -(tau_i . sigma_S n_S) = alpha_Qi tau_i.u_S.
-- Their alpha_Qi is the complete dimensional friction coefficient. Here it is
-- evaluated using the standard BJS scaling alpha_BJ*nu/sqrt(K_interface).
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
-- Qi et al. have no R_n (u.n)n term in eqs. (3)--(5) or weak form (7).
stokesBJS:set_normal_factor(0.0)

--------------------------------------------------------------------------------
-- Darcy equation through Membrane and Inter with a scaled pressure unknown
--------------------------------------------------------------------------------

darcyMembraneEq = ConvectionDiffusionFV1("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneDarcyMobility)

darcyInterEq = ConvectionDiffusionFV1("p_darcy_scaled", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterDarcyMobility)

-- Qi et al. impose zero Darcy flux on the exterior kidney boundary because
-- their porous domain couples an arterial tree to a venous tree. This tube has
-- only one lumen interface, so a fully sealed steady porous region could not
-- accept net leakage and would leave Darcy pressure without a reference. Keep
-- the validated OuterWall pressure sink; only the interface law is compared.
darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(
    outerPressureScaled, "p_darcy_scaled", "OuterWall"
)

coupledDisc = DomainDiscretization(approxSpace)
coupledDisc:add(stokesEq)
coupledDisc:add(stokesInlet)
coupledDisc:add(stokesOpenBoundary)
coupledDisc:add(stokesRim)
coupledDisc:add(stokesBJS)

-- Qi et al., eq. (4) and weak form (7):
--     -n_S . sigma_S n_S = p_D,
--     +(p_D, v_S.n_S)_Gamma.
-- This is normal-stress balance, not literal equality p_S=p_D.
-- NeumannBoundary takes the outward conservative momentum flux
-- -sigma_S n_S = p_D n_S.
local darcyPhysicalPressureValue =
    darcyPressureScale * darcyMembraneEq:value()
-- Keep every basis vector, product user data, and boundary discretization
-- alive for as long as coupledDisc exists.  Creating these only as loop-local
-- temporaries caused nondeterministic use-after-free crashes in
-- DataEvaluator::finish_elem_loop during matrix assembly.
local interfaceRows = {
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}
local interfaceTractionData = {}
local interfaceTractions = {}
for index, row in ipairs(interfaceRows) do
    interfaceTractionData[index] =
        interfaceTractionSign * darcyPhysicalPressureValue * row[2]
    local interfaceTraction = NeumannBoundary(row[1], "fv1")
    interfaceTraction:add(
        interfaceTractionData[index],
        "Apical", "Lumen"
    )
    interfaceTractions[index] = interfaceTraction
    coupledDisc:add(interfaceTraction)
end

local stokesVelocity = stokesEq:velocity()
darcyInterfaceFlux = NeumannBoundary("p_darcy_scaled", "fv1")
-- Evaluate the Stokes trace from the adjacent Membrane element, where the
-- velocity has a pinned extension but retains the shared Apical face DoFs.
-- NeumannBoundary converts the vector data to its outward normal component.
darcyInterfaceFlux:add(
    darcyEquationScale * interfaceFluxSign * stokesVelocity,
    "Apical", "Membrane"
)

-- Pin only non-physical extensions away from the constrained trace.
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(
    0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet"
)
inactiveStokes = DirichletBoundary()
for _, fct in ipairs({"u", "v", "w"}) do
    inactiveStokes:add(
        0.0, fct, "Membrane,Basolateral,MembraneWall"
    )
end
inactiveStokes:add(
    0.0, "p_stokes",
    "Membrane,Inter,Basolateral,OuterWall,MembraneWall"
)

coupledDisc:add(darcyMembraneEq)
coupledDisc:add(darcyInterEq)
coupledDisc:add(darcyInterfaceFlux)
coupledDisc:add(darcyOuterPressure)
coupledDisc:add(inactiveDarcy)
coupledDisc:add(inactiveStokes)
-- Uniform refinement of the hybrid prism/tetrahedron mesh creates a valid
-- nonconforming hierarchy at element transitions.  Replace every hanging P1
-- slave row by its one-sided interpolation relation on all GMG levels; without
-- this constraint r1 contains exact zero pressure rows at the inlet.
hangingNodeConstraints = OneSideP1Constraints()
coupledDisc:add(hangingNodeConstraints)
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- Qi et al. use right-preconditioned restarted GMRES(100). This comparison
-- uses GMRES(100) with geometric multigrid as the preconditioner, retaining
-- the scalable preconditioner required for the full nephron problem.

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    assert(solverSetup == "gmg", "Unknown solverSetup: " .. solverSetup)

    -- Split only the Stokes multiplier.  Velocity and Darcy pressure stay in
    -- the coupled elliptic forward block; a pressure mass matrix approximates
    -- the Stokes Schur complement.  OneSideP1Constraints is essential here so
    -- neither sliced block contains zero hanging-node rows.
    local gmgMainInverse = ILU()
    gmgMainInverse:set_damp(1.0)
    gmgMainInverse:set_inversion_eps(1.0e-30)
    local gmgPressureInverse = ILU()
    gmgPressureInverse:set_damp(1.0)
    gmgPressureInverse:set_inversion_eps(1.0e-30)
    local gmgSmoother = UzawaBase({"p_stokes"})
    gmgSmoother:set_forward_iter(gmgMainInverse)
    gmgSmoother:set_schur_iter(gmgPressureInverse)
    gmgSmoother:set_backward_iter(gmgMainInverse)

    local gmgSchurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
    gmgSchurMassEq:set_mass_scale(1.0)
    gmgSchurMassEq:set_diffusion(0.0)
    local gmgSchurMassDisc = DomainDiscretization(space)
    gmgSchurMassDisc:add(gmgSchurMassEq)
    gmgSchurMassDisc:add(OneSideP1Constraints())
    local gmgSchurMassOp = AssembledLinearOperator(gmgSchurMassDisc)
    gmgSmoother:set_schur_operator_update(
        gmgSchurMassOp, 1.0 / viscosity
    )

    return {
        type = "gmres",
        precond = {
            type = "gmg", -- alternative: "ilu" for a single-level solve
            approxSpace = space,
            discretization = disc,
            baseLevel = gmgBaseLevel,
            -- The outer/full-grid solve is iterative. SuperLU is restricted
            -- to the small level-0 saddle-point problem, where pivoting is
            -- needed for a reliable GMG coarse correction.
            baseSolver = "superlu",
            smoother = gmgSmoother,
            cycle = "V",
            preSmooth = 1,
            postSmooth = 1,
            rap = gmgUseRAP,
            gatheredBaseSolverIfAmbiguous = false,
            transfer = {
                type = "std",
                restrictionDamp = 1.0,
                prolongationDamp = 1.0,
                enableP1LagrangeOptimization = true
            }
        },
        convCheck = {
            type = "standard",
            iterations = 600,
            absolute = 1.0e-12,
            reduction = 1.0e-12,
            verbose = true
        }
    }
end

local function solveLinearSubproblem(space, disc, label)
    util.solver.defaults.approxSpace = space
    local matrix = AssembledLinearOperator(disc)
    local rhs = GridFunction(space)
    local solution = GridFunction(space)
    solution:set(0.0)
    disc:adjust_solution(solution)
    disc:assemble_linear(matrix, rhs)
    local componentNames = {
        "u", "v", "w", "p_stokes", "p_darcy_scaled"
    }
    local function trueComponentDefects()
        solution:set_consistent_storage_type()
        local defect = rhs:clone()
        VecAssign(defect, rhs)
        matrix:apply_sub(defect, solution)
        local norms = {}
        for selected = 1, #componentNames do
            local factors = {0.0, 0.0, 0.0, 0.0, 0.0}
            factors[selected] = 1.0
            local componentDefect = defect:clone()
            ScaleGF(componentDefect, defect, factors)
            norms[selected] = VecNorm(componentDefect)
        end
        return norms
    end
    local initialDefects = trueComponentDefects()
    local initialGlobalScaleSquared = 0.0
    for component = 1, #componentNames do
        initialGlobalScaleSquared = initialGlobalScaleSquared
            + initialDefects[component] * initialDefects[component]
    end
    local initialGlobalScale = math.sqrt(initialGlobalScaleSquared)
    local verifiedAbsolute = 1.0e-11
    local solverDesc = outerLinearSolverDescription(space, disc)
    -- The production path is always Krylov/GMG.  SuperLU may be used only by
    -- GMG on the tiny level-0 coarse problem, never on the refined matrix.
    local solver = GMRES(100)
    solver:set_preconditioner(
        util.solver.CreatePreconditioner(solverDesc.precond)
    )
    solver:set_convergence_check(
        util.solver.CreateConvCheck(solverDesc.convCheck)
    )
    print("=== " .. label .. " OUTER GMRES(100)/"
        .. solverSetup .. " convergence ===")
    solver:init(matrix, solution)
    local converged = false
    local solveCycles = maxKrylovCycles
    for restart = 1, solveCycles do
        print(string.format(
            "=== Solver cycle %d/%d ===", restart, solveCycles
        ))
        solver:apply(solution, rhs)
        local currentDefects = trueComponentDefects()
        converged = true
        for component = 1, #componentNames do
            -- Darcy has no independent right-hand side: its initial component
            -- defect is exactly zero and is generated only after Stokes--Darcy
            -- coupling becomes active. Use the global initial defect as the
            -- reference floor so this component receives the same meaningful
            -- relative tolerance as every non-zero initial component.
            local reference = math.max(
                initialDefects[component], initialGlobalScale
            )
            local relative = currentDefects[component] / reference
            print(string.format(
                "TRUE DEFECT %-8s abs=%.6e relative-to-reference=%.6e",
                componentNames[component], currentDefects[component], relative
            ))
            local target = verifiedAbsolute
            if currentDefects[component] ~= currentDefects[component]
                or currentDefects[component] > target then
                converged = false
            end
        end
        if converged then break end
    end
    assert(converged, label .. " restarted iterative/" .. solverSetup
        .. " failed the verified component tolerances")
    print("=== " .. label .. " OUTER solve converged ===")
    solution:set_consistent_storage_type()
    return solution
end

coupledState = solveLinearSubproblem(
    approxSpace, coupledDisc, "Monolithic Stokes--Darcy"
)
coupledState:set_consistent_storage_type()

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- This is a steady model, so no time discretization is required. The coupled
-- solution is initialized to zero inside solveLinearSubproblem before assembly.

-- Convert only the stored Darcy unknown back to physical pressure. Using this
-- state for postprocessing also makes darcyEq:gradient() a physical gradient.
physicalState = coupledState:clone()
ScaleGF(physicalState, coupledState, {
    1.0, 1.0, 1.0, 1.0, darcyPressureScale
})
physicalState:set_consistent_storage_type()

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- Writes lumen, membrane, and interstitial fields and reports pressure, flux,
-- conservation, and outlet/leakage split diagnostics.

darcyMembraneVelocity = DarcyVelocityLinker()
darcyMembraneVelocity:set_permeability(membranePermeability)
darcyMembraneVelocity:set_viscosity(viscosity)
darcyMembraneVelocity:set_density(1.0)
darcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
darcyMembraneVelocity:set_pressure_gradient(darcyMembraneEq:gradient())

darcyInterVelocity = DarcyVelocityLinker()
darcyInterVelocity:set_permeability(interPermeability)
darcyInterVelocity:set_viscosity(viscosity)
darcyInterVelocity:set_density(1.0)
darcyInterVelocity:set_gravity(ConstUserVector(0.0))
darcyInterVelocity:set_pressure_gradient(darcyInterEq:gradient())

local pressureTraceJumpL2 = L2Error(
    physicalState, "p_stokes", physicalState, "p_darcy_scaled", 2,
    "Apical"
)
local stokesInletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Inlet", "Lumen"
)
local stokesOutletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Outlet", "Lumen"
)
local stokesApicalFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Apical", "Lumen"
)
local darcyApicalFlow = IntegralNormalComponentOnManifold(
    darcyMembraneVelocity, physicalState, "Apical", "Membrane"
)
local darcyOuterWallFlow = IntegralNormalComponentOnManifold(
    darcyInterVelocity, physicalState, "OuterWall", "Inter"
)
-- Do not call either pressure Integral overload here.  The direct component
-- overload throws "trial space missing" for this mixed-support space, while
-- GridFunctionNumberData on a boundary face segfaults in this UG4 build.
-- Boundary pressure means are computed safely from the written VTU instead.
-- NavierStokesInflow assembles the analytic profile, whereas the solved P1
-- trace is only its nodal visualization.  Use the analytic FV1 inlet flux for
-- the physical balance and report the raw P1 result separately.
local prescribedInletFlowFV1 = IntegralNormalComponentOnManifold(
    "inletVelocity", physicalState, "Inlet", "Lumen", 0.0, 1
)
print(string.format(
    "Apical pressure-trace jump L2 (diagnostic, not constrained): %.6e",
    pressureTraceJumpL2
))
print(string.format(
    "OuterWall Darcy pressure: %.6e Pa; Stokes boundary means deferred to VTU postprocessing",
    outerPressurePa
))
print(string.format(
    "Lumen outward flows: prescribed inlet=%.6e, raw P1 inlet=%.6e, Outlet=%.6e, Apical=%.6e mm^3/s",
    prescribedInletFlowFV1, stokesInletFlow, stokesOutletFlow,
    stokesApicalFlow
))
local lumenBalance = prescribedInletFlowFV1
    + stokesOutletFlow + stokesApicalFlow
local lumenRelativeBalance = math.abs(lumenBalance)
    / math.max(math.abs(prescribedInletFlowFV1), 1.0e-30)
print(string.format(
    "Lumen FV1 physical balance=%.6e mm^3/s, relative=%.6e",
    lumenBalance, lumenRelativeBalance
))
local imposedInflow = -prescribedInletFlowFV1
local conservativeLeakFlow = imposedInflow - stokesOutletFlow
print(string.format(
    "Conservative flow split: Outlet=%.3f%%, leakage=%.3f%% of inlet",
    100.0 * stokesOutletFlow / imposedInflow,
    100.0 * conservativeLeakFlow / imposedInflow
))
print(string.format(
    "Raw Apical surface integral=%.3f%% of inlet; FV1 balance offset=%.6e mm^3/s",
    100.0 * stokesApicalFlow / imposedInflow,
    lumenBalance
))
print(string.format(
    "Apical diagnostic: Stokes=%.6e, raw Darcy-gradient=%.6e, sum=%.6e mm^3/s",
    stokesApicalFlow, darcyApicalFlow,
    stokesApicalFlow + darcyApicalFlow
))
print(string.format(
    "Darcy raw-gradient boundary flows: Apical=%.6e, OuterWall=%.6e, balance=%.6e mm^3/s",
    darcyApicalFlow, darcyOuterWallFlow,
    darcyApicalFlow + darcyOuterWallFlow
))
print(string.format(
    "Darcy assembled-boundary check: imposed Apical=%.6e, reconstructed OuterWall=%.6e, offset=%.6e mm^3/s",
    -stokesApicalFlow, darcyOuterWallFlow,
    -stokesApicalFlow + darcyOuterWallFlow
))
local darcyExitRelativeMismatch = math.abs(
    -stokesApicalFlow + darcyOuterWallFlow
) / math.max(math.abs(stokesApicalFlow), 1.0e-30)
local physicalValidationPassed = lumenRelativeBalance <= 1.0e-2
    and darcyExitRelativeMismatch <= 1.0e-2
print(string.format(
    "PHYSICAL VALIDATION %s: lumen relative balance=%.6e, Darcy exit relative mismatch=%.6e, tolerance=1.000000e-02",
    physicalValidationPassed and "PASS" or "FAIL",
    lumenRelativeBalance, darcyExitRelativeMismatch
))

local function vectorComponents(vectorData)
    local components = {}
    for component = 0, 2 do
        components[component + 1] = UserVectorEntryAdapter3d()
        components[component + 1]:set_vector(vectorData, component)
    end
    return components
end
local darcyMembraneVelocityComponents = vectorComponents(darcyMembraneVelocity)
local darcyInterVelocityComponents = vectorComponents(darcyInterVelocity)

stokesOut = VTKOutput()
stokesOut:clear_selection()
stokesOut:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
stokesOut:select_nodal("u", "u")
stokesOut:select_nodal("v", "v")
stokesOut:select_nodal("w", "w")
stokesOut:select_nodal("p_stokes", "stokes_pressure_kinematic")
stokesOut:print_subsets(
    vtk_file_name .. "_lumen", physicalState, "Lumen", 0, 0.0
)

local function writeDarcyOutput(suffix, subset, velocity, components)
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal("p_darcy_scaled", "darcy_pressure_kinematic")
    output:select_element(velocity, "darcy_velocity_mm_s")
    output:select_element(components[1], "u")
    output:select_element(components[2], "v")
    output:select_element(components[3], "w")
    output:print_subsets(
        vtk_file_name .. suffix, physicalState, subset, 0, 0.0
    )
end
writeDarcyOutput(
    "_membrane", "Membrane", darcyMembraneVelocity,
    darcyMembraneVelocityComponents
)
writeDarcyOutput(
    "_inter", "Inter", darcyInterVelocity,
    darcyInterVelocityComponents
)

local solve_elapsed_time = GetClockS() - solve_start_time
local total_elapsed_time = GetClockS() - total_start_time
print(string.format(
    "Monolithic Stokes--Darcy solve/postprocessing: %.6f mins",
    solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:             %.6f mins",
    total_elapsed_time / 60.0
))
