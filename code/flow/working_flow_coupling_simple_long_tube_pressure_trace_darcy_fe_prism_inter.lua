--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Straight projected-mesh validation candidate. Exact pressure continuity is
-- imposed by eliminating the stored Darcy trace rather than adding a
-- zero-diagonal multiplier unknown. Test scope on this workstation is r0--r2.
-- Solves a steady monolithic Stokes--Darcy system with Beavers--Joseph--
-- Saffman interface coupling. Stokes velocity/pressure are defined in the
-- lumen, scaled Darcy pressure is defined in membrane/interstitium, and the
-- reported Darcy pressure and velocity are converted back to physical units.
-- A weighted algebraic trace constraint imposes literal p_S = p_D on Apical,
-- combines the two conservation rows consistently, and retains Darcy scaling.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolves paths, loads UG4 utilities, and requires the discretization,
-- projector, and coarse-grid solver plugins used below.
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
gridS = "straight_nephron_prism_inter"
gridName = "runs/" .. gridS .. ".ugx"
numPreRefs = 0
-- Mesh-validated straight lumen: z in [-0.42, 0.42] mm, centered at the
-- origin, with the Apical wall at r = 0.015 mm.
inletRadius = 1.5e-2               -- [mm]
inletCenterX = 0.0                 -- [mm]
inletCenterY = 0.0                 -- [mm]
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
membranePermeabilityM2 = util.GetParamNumber(
    "-KmemM2", 1.0e-15, "membrane Darcy permeability [m^2]"
)
interPermeabilityM2 = util.GetParamNumber(
    "-KinterM2", 1.0e-14, "interstitial Darcy permeability [m^2]"
)
bjsInterfacePermeabilityM2 = util.GetParamNumber(
    "-KbjsM2", 1.0e-12, "effective BJS interface permeability [m^2]"
)
bjsAlpha = util.GetParamNumber(
    "-bjsAlpha", 0.2, "dimensionless BJS coefficient"
)
density = 1000.0                  -- [kg/m^3]
outerPressurePa = util.GetParamNumber(
    "-outerPressurePa", 4.6, "interstitial OuterWall pressure [Pa]"
)
outletPressurePa = util.GetParamNumber(
    "-outletPressurePa", 0.0, "prescribed lumen outlet normal stress [Pa]"
)
darcyPressureScale = util.GetParamNumber(
    "-pressureScale", 1.0e4,
    "algebraic Darcy-pressure scale; output is rescaled physically"
)
normalVelocityStabilization = util.GetParamNumber(
    "-normalImpedance", 1.0e3, "Apical numerical normal impedance"
)
solverSetup = "gmg"
outerKrylovType = "bicgstab"      -- GMRES produced NaNs with this preconditioner
gmgSmootherType = "uzawa"          -- Stokes pressure is the saddle component
gmgUseRAP = true                  -- preserve the constrained operator by Galerkin products
interfaceTractionSign = 1.0
maxKrylovCycles = 8
krylovCycleReduction = 1.0e-4    -- per-cycle algebraic reduction target
gmgBaseLevel = 0                  -- retain the loaded coarse mesh as the GMG base
darcyDiscretization = util.GetParam(
    "-darcyDisc", "fe", "Darcy discretization: fe or fv1"
)

-- MPI runs create .pvtu master files and one .vtu file per rank.
local flowModelTag = bstokes and "stokes" or "navier_stokes"
vtk_file_name = util.GetParam(
    "-output",
    "Results_final/" .. gridS .. "_" .. flowModelTag
        .. "_darcy_bjs_pressure_trace_constraint_" .. solverSetup
        .. "_ref" .. numRefs,
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
assert(normalVelocityStabilization >= 0.0,
    "normalVelocityStabilization must be non-negative")
assert(darcyDiscretization == "fe" or darcyDiscretization == "fv1",
    "darcyDisc must be 'fe' or 'fv1'")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print(string.format(
    "Model: %s, vmax=%.6e mm/s, refinements=%d",
    bstokes and "Stokes" or "Navier--Stokes", vmax, numRefs
))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Loads the straight-tube mesh together with its serialized NeuriteProjector,
-- validates the projector data and subsets, then refines the projected
-- geometry before creating the coupled space.

local function file_contains(path, needle)
    local file = assert(io.open(path, "r"), "Cannot open grid file: " .. path)
    for line in file:lines() do
        if string.find(line, needle, 1, true) then
            file:close()
            return true
        end
    end
    file:close()
    return false
end

assert(file_contains(gridName, '<projection_handler'),
    "The selected UGX has no serialized ProjectionHandler: " .. gridName)
assert(file_contains(gridName, '<projector type="neurite"'),
    "The selected UGX has no serialized NeuriteProjector: " .. gridName)
assert(file_contains(gridName, '<vertex_attachment name="npSurfParams"'),
    "The selected UGX has no npSurfParams refinement data: " .. gridName)

print("Projected straight-tube grid: " .. gridName)

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Membrane", "Inter", "Apical", "Basolateral",
     "OuterWall", "MembraneWall", "Inlet", "Outlet", "TripleJunction"}),
    "The domain is missing a required Stokes--Darcy subset")

-- Distributed construction of this coupled constrained GMG hierarchy produces
-- either MPI_ERR_TRUNCATE at level 0 or NaNs when distribution starts on an
-- intermediate level.  Refine with the serialized projector first, then
-- distribute the completed hierarchy at its finest level.
balancer.firstDistLvl = numRefs
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
    membraneDarcyMobility * darcyPressureScale
local scaledInterDarcyMobility = interDarcyMobility * darcyPressureScale
local bjsFriction =
    bjsAlpha * viscosity / math.sqrt(bjsInterfacePermeability)
local outerPressurePhysical = outerPressurePa / (density * 1.0e-6)
local outerPressureScaled = outerPressurePhysical / darcyPressureScale
local outletPressureKinematic = outletPressurePa / (density * 1.0e-6)

print(string.format(
    "K_mem=%.6e mm^2, K_inter=%.6e mm^2, pressure scale=%.6e",
    membranePermeability, interPermeability, darcyPressureScale
))
print(string.format(
    "Scaled Darcy diffusion: membrane=%.6e, inter=%.6e; BJS friction=%.6e; normal impedance=%.6e; OuterWall pressure=%.6e Pa",
    scaledMembraneDarcyMobility, scaledInterDarcyMobility,
    bjsFriction, normalVelocityStabilization, outerPressurePa
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

-- BJS--Saffman: (sigma_S n)_t = -alpha*nu/sqrt(K) u_S,t.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
-- A small positive normal impedance damps interface-grid modes. This is far
-- below the 5e6 membrane resistance in the validated baseline and is being
-- tested explicitly for GMG robustness.
stokesBJS:set_normal_factor(normalVelocityStabilization)

--------------------------------------------------------------------------------
-- Darcy equation through Membrane and Inter with a scaled pressure unknown
--------------------------------------------------------------------------------

local DarcyDiscretization = darcyDiscretization == "fv1"
    and ConvectionDiffusionFV1 or ConvectionDiffusionFE
darcyMembraneEq = DarcyDiscretization("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneDarcyMobility)

darcyInterEq = DarcyDiscretization("p_darcy_scaled", "Inter")
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

-- Laptev (2003), eq. (1.152), and Jaeger--Mikelic--Neuss give the Darcy
-- pressure as the free-fluid normal stress (plus a microstructure-dependent
-- boundary-layer shear correction).  For a symmetric, isotropic interface
-- that correction is zero.  Thus sigma_S n_S = -p_D n_S; this is not literal
-- equality of the two pressure unknowns.  NeumannBoundary takes the outward
-- conservative momentum flux -sigma_S n_S = p_D n_S.
local darcyPhysicalPressureValue =
    darcyPressureScale * darcyMembraneEq:value()
for _, row in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    local interfaceTraction = NeumannBoundary(row[1], "fv1")
    interfaceTraction:add(
        interfaceTractionSign * darcyPhysicalPressureValue * row[2],
        "Apical", "Lumen"
    )
    coupledDisc:add(interfaceTraction)
    if outletPressurePa ~= 0.0 then
        local outletPressureTraction = NeumannBoundary(row[1], "fv1")
        outletPressureTraction:add(
            outletPressureKinematic * row[2], "Outlet", "Lumen"
        )
        coupledDisc:add(outletPressureTraction)
    end
end

-- Eliminate p_D,stored on Apical by p_D,stored=p_S/scale. The eliminated
-- Darcy balance is added to the Stokes continuity row without row scaling;
-- only the eliminated Darcy column carries 1/scale.
pressureContinuity = PressureTraceConstraint(
    "p_stokes", "p_darcy_scaled", "Apical,TripleJunction"
)
pressureContinuity:set_pressure_scale(darcyPressureScale)

-- The row merge above is the normal-flux coupling: Stokes continuity plus
-- Darcy conservation is one interface balance. Do not also prescribe u_S as
-- a Darcy Neumann flux, which would count the Stokes normal flux twice.
local stokesVelocity = stokesEq:velocity()

-- Coupled user data must see both fields on the elements adjacent to Apical.
-- Pin only their non-physical extensions; never pin the Apical trace itself.
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(
    0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet"
)
inactiveStokes = DirichletBoundary()
inactiveStokes:add(
    0.0, "p_stokes", "Membrane,Basolateral,MembraneWall"
)

coupledDisc:add(darcyMembraneEq)
coupledDisc:add(darcyInterEq)
coupledDisc:add(pressureContinuity)
coupledDisc:add(darcyOuterPressure)
coupledDisc:add(inactiveDarcy)
coupledDisc:add(inactiveStokes)
-- The hybrid prism/tetrahedron hierarchy is nonconforming at transition
-- faces. Replace hanging P1 slave rows by their interpolation relations on
-- every level used by the pressure trace constraint and GMG.
hangingNodeConstraints = OneSideP1Constraints()
coupledDisc:add(hangingNodeConstraints)
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- Solves the monolithic system with restarted BiCGStab and geometric
-- multigrid while printing algebraic and componentwise true defects.

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    assert(solverSetup == "gmg", "Unknown solverSetup: " .. solverSetup)

    -- Only p_stokes is a saddle-point component. Darcy pressure remains in the
    -- main elliptic block together with velocity and the constrained trace.
    assert(gmgSmootherType == "uzawa",
        "Unknown gmgSmootherType: " .. tostring(gmgSmootherType))
    -- Line smoothing is needed for the strongly anisotropic Darcy modes on
    -- the refined straight-tube hierarchy. Point ILU stagnates at r2.
    local gmgMainInverse = LineGaussSeidel(space)
    gmgMainInverse:set_num_steps(1, 1, 1, 1, 1, 1)
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
    gmgSmoother:set_schur_operator_update(gmgSchurMassOp, 1.0 / viscosity)

    return {
        type = outerKrylovType,
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
            preSmooth = 2,
            postSmooth = 2,
            rap = gmgUseRAP,
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
            absolute = 1e-14,
            reduction = krylovCycleReduction,
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
    local verifiedReduction = 1.0e-8
    local verifiedAbsolute = 1.0e-12
    local solverDesc = outerLinearSolverDescription(space, disc)
    local solver = util.solver.CreateLinearSolver(solverDesc)
    print("=== " .. label .. " OUTER " .. outerKrylovType .. "/"
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
            local target = math.max(
                verifiedAbsolute,
                reference * verifiedReduction
            )
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
local inletPressureMean = Integral(
    physicalState, "p_stokes", "Inlet"
) / Integral(1.0, physicalState, "Inlet")
local outletPressureMean = Integral(
    physicalState, "p_stokes", "Outlet"
) / Integral(1.0, physicalState, "Outlet")
print(string.format(
    "Apical pressure-trace jump L2 (constraint includes TripleJunction rings): %.6e",
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
    "Monolithic Stokes--Darcy %s runtime: %.6f mins",
    solverSetup, solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:             %.6f mins",
    total_elapsed_time / 60.0
))
