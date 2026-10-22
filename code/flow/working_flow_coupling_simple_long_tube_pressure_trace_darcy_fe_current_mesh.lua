--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- This script solves one steady, monolithic Stokes--Darcy problem:
--
--   * Stokes flow is assembled only in the Lumen volume.
--   * Darcy flow is assembled in the Membrane and Inter volumes.
--   * The Apical surface couples both models through pressure continuity,
--     normal mass conservation, and the Beavers--Joseph--Saffman tangential
--     slip condition.
--   * p_stokes is the physical kinematic pressure. p_darcy_scaled stores
--     p_D / darcyPressureScale for matrix conditioning; output multiplies it
--     back by darcyPressureScale. The scaling does not change the physics.
--
-- PressureTraceConstraint eliminates the duplicate Darcy trace on Apical and
-- imposes p_stokes = darcyPressureScale * p_darcy_scaled. It also merges the
-- Stokes and Darcy conservation rows, which enforces normal-flux continuity
-- without introducing a zero-diagonal Lagrange-multiplier field.
--
-- This is the preserved formulation that solved numRefs=2 and produced the
-- pressure_trace_darcy_fe_ref2 output. The known refined-grid Darcy-pressure
-- oscillation is intentionally not hidden or filtered in this script.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolve the Model directory from this file instead of relying on the shell's
-- current directory. Load only the UG4 utility scripts and plugins used by the
-- mesh hierarchy, finite-volume Stokes operator, finite-element Darcy operator,
-- NeuriteProjector, and exact GMG coarse solve.
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
-- units, boundary data, algebraic scaling, and GMG tolerances live here so a
-- run can be reproduced from this file alone.

dim = 3

numRefs = util.GetParamNumber("-numRefs", 2, "number of uniform refinements")
vmax = 3.0                          -- maximum inlet velocity [mm/s]

-- Geometry and physical parameters: edit values here, not on the command line.
gridS = "stright_tube"
gridName = "runs/" .. gridS .. ".ugx"
-- Physical lumen radius.  The ring at r = 5e-3 mm is the inner O-grid
-- transition, not the wall; the lumen/TripleJunction rim is at r = 1e-2 mm.
inletRadius = 1.0e-2               -- [mm]
inletCenterX = 2.286625           -- [mm]
inletCenterY = 3.331125           -- [mm]
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
membranePermeabilityM2 = 1.0e-15  -- membrane Darcy permeability [m^2]
interPermeabilityM2 = 1.0e-14     -- interstitial Darcy permeability [m^2]
bjsInterfacePermeabilityM2 = 1.0e-12 -- effective BJS slip permeability [m^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- [kg/m^3]
outerPressurePa = 4.6             -- [Pa], calibrated interstitial back-pressure
darcyPressureScale = 1.0e4        -- algebraic only; output is rescaled physically
normalVelocityStabilization = 1.0e3
                                    -- Apical normal impedance [mm/s]
gmgBaseLevel = 0                  -- retain the loaded coarse mesh as the GMG base

-- MPI runs create .pvtu master files and one .vtu file per rank.
vtk_file_name = util.GetParam(
    "-output",
    "Results_final/" .. gridS
        .. "_stokes_darcy_bjs_pressure_trace_constraint_gmg"
        .. "_ref" .. numRefs,
    "VTK output prefix"
)


InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- The UGX contains the coarse mesh, subset layout, NeuriteProjector, and the
-- npSurfParams attachment required to place new surface vertices.

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Membrane", "Inter", "Apical", "Basolateral",
     "OuterWall", "MembraneWall", "Inlet", "Outlet", "TripleJunction"}),
    "The domain is missing a required Stokes--Darcy subset")

-- UG4 commit 24ff76c4c broadcasts the serialized ProjectionHandler when the
-- grid is loaded on rank 0. Distribute the base grid before refinement so MPI
-- ranks do not retain the full refined mesh as ghosts.
balancer.firstDistLvl = 0
-- Keep all load-balancing choices in this script; do not expose UG4's
-- additional balancing command-line parameters.
balancer.parametersParsed = true
balancer.RefineAndRebalanceDomain(dom, numRefs)

--------------------------------------------------------------------------------
-- Approximation-space fields
--------------------------------------------------------------------------------
-- All five unknowns belong to one ApproximationSpace and therefore one global
-- matrix. Their subset support is deliberately wider than their physical PDE
-- domains near Apical: UG4 user-data evaluation needs both pressure traces on
-- elements adjacent to the interface. The non-physical extensions are pinned
-- to zero below and do not acquire Stokes or Darcy volume equations.

approxSpace = ApproximationSpace(dom)
local allFlowSubsets =
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet,TripleJunction"
local stokesVelocitySubsetNames = {
    "Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"
}
local stokesPressureAndInterfaceSubsetNames = {
    "Lumen", "Membrane", "Apical", "Basolateral", "MembraneWall",
    "Inlet", "Outlet", "TripleJunction"
}
approxSpace:add_fct(
    {"u", "v", "w"}, "Lagrange", 1, stokesVelocitySubsetNames
)
-- p_stokes has a pinned membrane-side extension because the inner-boundary
-- pressure constraint is evaluated from both volumes adjacent to Apical. No
-- Stokes momentum or continuity operator is assembled in Membrane.
approxSpace:add_fct(
    {"p_stokes"}, "Lagrange", 1, stokesPressureAndInterfaceSubsetNames
)
approxSpace:add_fct(
    "p_darcy_scaled", "Lagrange", 1, allFlowSubsets
)
approxSpace:init_levels()
approxSpace:init_top_surface()

-- Print UG4's global distributed-domain statistics followed by the algebraic
-- DoF statistics. The first table reports the mesh at every refinement level;
-- the second reports the unknown count contributed by each solution field.
print("=== MESH STATISTICS ===")
print(string.format("Mesh file: %s", gridName))
print(string.format("Requested global refinements: %d", numRefs))
print(string.format("MPI processes: %d", NumProcs()))
print(dom:domain_info():to_string())
print("=== APPROXIMATION-SPACE / DOF STATISTICS ===")
approxSpace:print_statistic("component")

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Convert SI permeability to mm^2, construct the parabolic inlet profile, and
-- assemble the Stokes, Darcy, and interface equations. Pressures inside UG4 are
-- kinematic pressures [mm^2/s^2]; multiplication by density*1e-6 converts them
-- to Pa for reporting.

local membranePermeability = membranePermeabilityM2 * 1.0e6 -- [mm^2]
local interPermeability = interPermeabilityM2 * 1.0e6       -- [mm^2]
local bjsInterfacePermeability = bjsInterfacePermeabilityM2 * 1.0e6
local membraneDarcyMobility = membranePermeability / viscosity
local interDarcyMobility = interPermeability / viscosity
local scaledMembraneDarcyMobility =
    membraneDarcyMobility * darcyPressureScale
local scaledInterDarcyMobility = interDarcyMobility * darcyPressureScale
local bjsFriction =
    bjsAlpha * viscosity / math.sqrt(bjsInterfacePermeability)
local outerPressurePhysical = outerPressurePa / (density * 1.0e-6)
local outerPressureScaled = outerPressurePhysical / darcyPressureScale

function inletVelocity(x, y, z, t)
    -- The inlet is the z-normal circular cap of this reference mesh. Clamp the
    -- rim to zero so the prescribed profile is compatible with the wall/rim BC.
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

-- Stokes equation in Lumen
--------------------------------------------------------------------------------
-- FV1 uses equal-order P1 velocity and pressure. The fields/cor stabilization
-- suppresses the equal-order pressure null mode. set_laplace(false) selects the
-- symmetric-gradient stress tensor needed by the BJS traction condition.

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
stokesEq:set_stokes(true)
stokesEq:set_laplace(false) -- symmetric-gradient stress required for BJS
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_stabilization("fields", "cor")

-- Prescribe the full velocity at Inlet. Outlet receives the standard zero
-- normal-stress condition. The same open-boundary discretization is required
-- on Apical to assemble the boundary-normal flux in Stokes continuity.
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

-- BJS--Saffman tangential law:
--     (sigma_S n)_t = -alpha*nu/sqrt(K_interface) * u_S,t.
-- The normal factor is a Robin impedance used by the validated formulation to
-- damp an otherwise weak interface-normal mode; it does not replace the normal
-- conservation equation imposed by PressureTraceConstraint.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalVelocityStabilization)

--------------------------------------------------------------------------------
-- Darcy equation through Membrane and Inter with a scaled pressure unknown
--------------------------------------------------------------------------------
-- With q_D = -(K/nu) grad(p_D), substituting
-- p_D = darcyPressureScale*p_darcy_scaled multiplies the diffusion coefficient
-- by the same scale. Separate coefficients represent the lower-permeability
-- membrane and higher-permeability interstitium. A physical pressure is fixed
-- only on OuterWall; Basolateral is an internal Membrane/Inter interface.

darcyMembraneEq = ConvectionDiffusionFE("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneDarcyMobility)

darcyInterEq = ConvectionDiffusionFE("p_darcy_scaled", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterDarcyMobility)

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

-- Normal traction balance supplies Darcy pressure to the Stokes momentum rows.
-- The neglected viscous normal-stress correction leaves sigma_S*n_S=-p_D*n_S.
-- NeumannBoundary uses the outward conservative momentum flux -sigma_S*n_S,
-- so the applied vector traction is +p_D times each Cartesian basis vector.
local darcyPhysicalPressureValue =
    darcyPressureScale * darcyMembraneEq:value()
for _, row in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    local interfaceTraction = NeumannBoundary(row[1], "fv1")
    interfaceTraction:add(
        darcyPhysicalPressureValue * row[2],
        "Apical", "Lumen"
    )
    coupledDisc:add(interfaceTraction)
end

-- Eliminate the duplicate Darcy trace on Apical:
--     p_darcy_scaled = p_stokes / darcyPressureScale.
-- The eliminated Darcy conservation row is added to the Stokes continuity row.
-- Consequently this one operation imposes physical pressure continuity and
-- balances the Stokes and Darcy normal fluxes without double-counting either.
pressureContinuity = PressureTraceConstraint(
    "p_stokes", "p_darcy_scaled", "Apical"
)
pressureContinuity:set_pressure_scale(darcyPressureScale)

-- This velocity handle is used only for flux diagnostics and VTK output. No
-- additional Darcy Neumann condition is added on Apical because the merged
-- conservation row already supplies the normal-flux coupling.
local stokesVelocity = stokesEq:velocity()

-- Coupled user data must see both fields on the elements adjacent to Apical.
-- Pin only their non-physical extensions; never pin the Apical trace itself.
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(
    0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet"
)
inactiveStokes = DirichletBoundary()
inactiveStokes:add(
    0.0, "p_stokes", "Membrane,Basolateral,MembraneWall,TripleJunction"
)

coupledDisc:add(darcyMembraneEq)
coupledDisc:add(darcyInterEq)
coupledDisc:add(pressureContinuity)
coupledDisc:add(darcyOuterPressure)
coupledDisc:add(inactiveDarcy)
coupledDisc:add(inactiveStokes)
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- The coupled matrix is indefinite because p_stokes is the Stokes
-- incompressibility multiplier. BiCGStab is therefore preconditioned by GMG
-- with an Uzawa smoother. Uzawa treats p_stokes as the Schur component while
-- velocity and the elliptic Darcy pressure stay in the main block.
--
-- Each V-cycle uses two pre/post smoothing passes, RAP/Galerkin coarse
-- operators (important for the algebraically constrained Apical rows), and a
-- pivoted SuperLU solve only on the small level-0 problem. The large finest
-- system is always solved iteratively.

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    -- ILU approximates the main-block and pressure-mass inverses used by the
    -- forward, Schur, and backward stages of one symmetric Uzawa sweep.
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

    -- The inverse pressure mass matrix approximates the Stokes Schur
    -- complement. Scaling by 1/nu gives the correct viscous Stokes magnitude.
    local gmgSchurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
    gmgSchurMassEq:set_mass_scale(1.0)
    gmgSchurMassEq:set_diffusion(0.0)
    local gmgSchurMassDisc = DomainDiscretization(space)
    gmgSchurMassDisc:add(gmgSchurMassEq)
    local gmgSchurMassOp = AssembledLinearOperator(gmgSchurMassDisc)
    gmgSmoother:set_schur_operator_update(gmgSchurMassOp, 1.0 / viscosity)

    return {
        type = "bicgstab",
        precond = {
            type = "gmg",
            approxSpace = space,
            discretization = disc,
            baseLevel = gmgBaseLevel,
            -- The outer/full-grid solve is iterative. SuperLU is restricted
            -- to the small level-0 saddle-point problem, where pivoting is
            -- needed for a reliable GMG coarse correction.
            baseSolver = "superlu",
            smoother = gmgSmoother,
            cycle = "V",
            preSmooth = 2,
            postSmooth = 2,
            rap = true,
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
            -- Global BiCGStab convergence criteria are literal here so the
            -- complete stopping rule is visible in this solver block.
            iterations = 400,
            absolute = 1.0e-8,
            reduction = 1.0e-8,
            verbose = true
        }
    }
end

local function solveLinearSubproblem(space, disc, label)
    -- Assemble once because the present Stokes--Darcy problem is linear and
    -- steady. The solution starts at zero after essential BC adjustment.
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
    -- Recompute r=b-Ax from the assembled operator. This guards against a
    -- recursive Krylov residual appearing converged while one physical block
    -- remains inaccurate.
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
    local solverDesc = outerLinearSolverDescription(space, disc)
    local solver = util.solver.CreateLinearSolver(solverDesc)
    print("=== " .. label .. " OUTER bicgstab/gmg convergence ===")
    print("=== CONVERGENCE CRITERIA ===")
    print("BiCGStab maximum iterations per cycle: 400")
    print("BiCGStab absolute global defect:       <= 1.000000e-08")
    print("BiCGStab relative global reduction:    <= 1.000000e-08")
    print("Maximum restarted solver cycles:       8")
    print("Each recomputed component defect:      <= 1.000000e-08")
    print("A cycle passes only when every component is PASS.")
    print("=== INITIAL TRUE COMPONENT DEFECTS ===")
    for component = 1, #componentNames do
        print(string.format(
            "INITIAL DEFECT %-14s abs=%.6e",
            componentNames[component], initialDefects[component]
        ))
    end
    print(string.format(
        "INITIAL COMBINED DEFECT              abs=%.6e",
        initialGlobalScale
    ))
    solver:init(matrix, solution)
    local converged = false
    for restart = 1, 8 do
        print(string.format(
            "=== Solver cycle %d/8 ===", restart
        ))
        solver:apply(solution, rhs)
        local currentDefects = trueComponentDefects()
        local currentGlobalSquared = 0.0
        for component = 1, #componentNames do
            currentGlobalSquared = currentGlobalSquared
                + currentDefects[component] * currentDefects[component]
        end
        local currentGlobal = math.sqrt(currentGlobalSquared)
        local globalRelative = currentGlobal / initialGlobalScale
        print(string.format(
            "TRUE COMBINED DEFECT initial=%.6e current=%.6e "
                .. "relative=%.6e absolute_target=1.000000e-08 "
                .. "relative_target=1.000000e-08",
            initialGlobalScale, currentGlobal, globalRelative
        ))
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
            local passed = currentDefects[component] <= 1.0e-8
            print(string.format(
                "TRUE DEFECT %-14s initial=%.6e current=%.6e "
                    .. "relative=%.6e target=1.000000e-08 %s",
                componentNames[component], initialDefects[component],
                currentDefects[component], relative,
                passed and "PASS" or "FAIL"
            ))
            if currentDefects[component] ~= currentDefects[component]
                or not passed then
                converged = false
            end
        end
        if converged then break end
    end
    assert(converged, label .. " restarted iterative/gmg"
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
-- No time operator is assembled because the requested model is steady. The
-- zero algebraic initial guess is created inside solveLinearSubproblem. This
-- block converts the stored Darcy unknown to physical kinematic pressure for
-- all diagnostics and VTK fields; Stokes velocity and pressure are unchanged.

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
-- DarcyVelocityLinker evaluates q_D=-(K/nu)grad(p_D) separately with membrane
-- and interstitial permeability. The three VTK datasets deliberately contain
-- physical pressure and velocity values, not the internally scaled unknown.
-- Integral diagnostics report pressure continuity, global lumen mass balance,
-- and the fraction leaving through Outlet versus Apical leakage.

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
local inletPressureMean = Integral(
    physicalState, "p_stokes", "Inlet"
) / Integral(1.0, physicalState, "Inlet")
local outletPressureMean = Integral(
    physicalState, "p_stokes", "Outlet"
) / Integral(1.0, physicalState, "Outlet")
print(string.format(
    "Apical pressure-trace jump L2 (constraint excludes TripleJunction rings): %.6e",
    pressureTraceJumpL2
))
print(string.format(
    "Mean Stokes pressure: Inlet=%.6e Pa, Outlet=%.6e Pa; OuterWall=%.6e Pa",
    density * 1.0e-6 * inletPressureMean,
    density * 1.0e-6 * outletPressureMean, outerPressurePa
))
print(string.format(
    "Lumen outward flows: Inlet=%.6e, Outlet=%.6e, Apical=%.6e, balance=%.6e mm^3/s",
    stokesInletFlow, stokesOutletFlow, stokesApicalFlow,
    stokesInletFlow + stokesOutletFlow + stokesApicalFlow
))
local imposedInflow = -stokesInletFlow
local conservativeLeakFlow = imposedInflow - stokesOutletFlow
print(string.format(
    "Conservative flow split: Outlet=%.3f%%, leakage=%.3f%% of inlet",
    100.0 * stokesOutletFlow / imposedInflow,
    100.0 * conservativeLeakFlow / imposedInflow
))
print(string.format(
    "Raw Apical surface integral=%.3f%% of inlet; quadrature/discrete balance offset=%.6e mm^3/s",
    100.0 * stokesApicalFlow / imposedInflow,
    stokesInletFlow + stokesOutletFlow + stokesApicalFlow
))
print(string.format(
    "Apical diagnostic: Stokes=%.6e, raw Darcy-gradient=%.6e, sum=%.6e mm^3/s",
    stokesApicalFlow, darcyApicalFlow,
    stokesApicalFlow + darcyApicalFlow
))

local function vectorComponents(vectorData)
    -- VTK receives both the vector and named scalar components so ParaView can
    -- display arrows or inspect u, v, and w independently.
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
    -- Darcy velocity is element data because it is reconstructed from the P1
    -- pressure gradient, which is elementwise constant for affine elements.
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
    "Monolithic Stokes--Darcy GMG solve and postprocessing: %.6f mins",
    solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:                              %.6f mins",
    total_elapsed_time / 60.0
))
