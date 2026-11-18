--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- EXPERIMENTAL: uses one physical pressure unknown in Lumen, Apical,
-- Membrane, and Inter, so literal pressure continuity is built into the space.
-- Solves a steady monolithic Stokes--Darcy system with Beavers--Joseph--
-- Saffman interface coupling. Navier--Stokes is assembled only in Lumen and
-- Darcy only in Membrane/Inter. Their pressure rows meet on Apical without a
-- multiplier, trace elimination, or duplicated interface flux equation.

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
gridS = "stright_tube"
gridName = "runs/" .. gridS .. ".ugx"
numPreRefs = 0
-- Physical lumen radius.  The ring at r = 5e-3 mm is the inner O-grid
-- transition, not the wall; the lumen/MembraneWall rim is at r = 1e-2 mm.
inletRadius = 1.0e-2               -- [mm]
inletCenterX = 2.286625           -- [mm]
inletCenterY = 3.331125           -- [mm]
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
membranePermeabilityM2 = 1.0e-15  -- membrane Darcy permeability [m^2]
interPermeabilityM2 = 1.0e-14     -- interstitial Darcy permeability [m^2]
bjsInterfacePermeabilityM2 = 1.0e-12 -- effective BJS slip permeability [m^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- [kg/m^3]
outerPressurePa = 0.0             -- [Pa], no imposed interstitial back-pressure
globalPressureScale = 1.0          -- physical shared pressure; no domainwise scaling
solverSetup = "gmg"
outerKrylovType = "bicgstab"      -- GMRES produced NaNs with this preconditioner
gmgSmootherType = "uzawa"          -- Stokes pressure is the saddle component
gmgUseRAP = true                  -- preserve the constrained operator by Galerkin products
interfaceTractionSign = 1.0
maxKrylovCycles = 8
gmgBaseLevel = 0                  -- retain the loaded coarse mesh as the GMG base

-- MPI runs create .pvtu master files and one .vtu file per rank.
local flowModelTag = bstokes and "stokes" or "navier_stokes"
vtk_file_name = util.GetParam(
    "-output",
    "Results_final/" .. gridS .. "_" .. flowModelTag
        .. "_darcy_bjs_shared_pressure_" .. solverSetup
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

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print(string.format(
    "Model: %s, vmax=%.6e mm/s, refinements=%d",
    bstokes and "Stokes" or "Navier--Stokes", vmax, numRefs
))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Loads the straight tube and its embedded refinement projectors, validates
-- its subsets, refines the domain, and creates one coupled space.

local function file_contains(path, needle)
    -- Confirms that refinement-critical projector metadata exists before load.
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
    "The straight-tube grid has no serialized NeuriteProjector: " .. gridName)
assert(file_contains(gridName, '<vertex_attachment name="npSurfParams"'),
    "The straight-tube grid has no npSurfParams refinement data: " .. gridName)

print("Neurite-projector straight-tube grid: " .. gridName)

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Membrane", "Inter", "Apical", "Basolateral",
     "OuterWall", "MembraneWall", "Inlet", "Outlet"}),
    "The domain is missing a required Stokes--Darcy subset")

-- UG4 commit 24ff76c4c broadcasts the serialized ProjectionHandler when the
-- grid is loaded on rank 0, so the former double-load registration workaround
-- is no longer required. Distribution is currently delayed until the finest
-- level because distributed construction of this coupled GMG hierarchy still
-- makes the outer Krylov iteration stagnate; this is separate from projector
-- deserialization and is checked by the true-defect verifier below.
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
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet"
local stokesVelocitySubsetNames = {
    "Lumen", "Apical", "Inlet", "Outlet", "MembraneWall"
}
approxSpace:add_fct(
    {"u", "v", "w"}, "Lagrange", 1, stokesVelocitySubsetNames
)
approxSpace:add_fct(
    "p_pressure", "Lagrange", 1, allFlowSubsets
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
local scaledMembraneDarcyMobility = membraneDarcyMobility
local scaledInterDarcyMobility = interDarcyMobility
local bjsFriction =
    bjsAlpha * viscosity / math.sqrt(bjsInterfacePermeability)
local outerPressurePhysical = outerPressurePa / (density * 1.0e-6)
local outerPressureScaled = outerPressurePhysical

print(string.format(
    "K_mem=%.6e mm^2, K_inter=%.6e mm^2, pressure scale=%.6e",
    membranePermeability, interPermeability, globalPressureScale
))
print(string.format(
    "Scaled Darcy diffusion: membrane=%.6e, inter=%.6e; BJS friction=%.6e; OuterWall pressure=%.6e Pa",
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

stokesEq = NavierStokesFV1({"u", "v", "w", "p_pressure"}, {"Lumen"})
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
stokesRim:add("MembraneWall")

-- BJS--Saffman: (sigma_S n)_t = -alpha*nu/sqrt(K) u_S,t.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
-- Literal p_S=p_D with the viscous normal-stress correction neglected leaves
-- no additional normal Robin term; retaining it would impose R_n*u_n=0.
stokesBJS:set_normal_factor(0.0)

--------------------------------------------------------------------------------
-- Darcy equation through Membrane and Inter with the shared physical pressure
--------------------------------------------------------------------------------

darcyMembraneEq = ConvectionDiffusionFV1("p_pressure", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneDarcyMobility)

darcyInterEq = ConvectionDiffusionFV1("p_pressure", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterDarcyMobility)

darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(
    outerPressureScaled, "p_pressure", "OuterWall"
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
local darcyPhysicalPressureValue = darcyMembraneEq:value()
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
end

-- The single pressure row on Apical receives the Stokes-continuity and Darcy
-- conservation contributions. Their sum imposes normal-flux conservation;
-- adding a separate Darcy Neumann term here would impose it twice.
local stokesVelocity = stokesEq:velocity()

coupledDisc:add(darcyMembraneEq)
coupledDisc:add(darcyInterEq)
coupledDisc:add(darcyOuterPressure)
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- Solves the monolithic system with restarted BiCGStab and geometric
-- multigrid while printing algebraic and componentwise true defects.

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    assert(solverSetup == "gmg", "Unknown solverSetup: " .. solverSetup)

    -- The shared pressure is the only Schur component. Its approximate Schur
    -- operator contains both the lumen pressure mass and porous diffusion.
    assert(gmgSmootherType == "uzawa",
        "Unknown gmgSmootherType: " .. tostring(gmgSmootherType))
    local gmgMainInverse = ILU()
    gmgMainInverse:set_damp(1.0)
    gmgMainInverse:set_inversion_eps(1.0e-30)
    local gmgPressureInverse = ILU()
    gmgPressureInverse:set_damp(1.0)
    gmgPressureInverse:set_inversion_eps(1.0e-30)
    local gmgSmoother = UzawaBase({"p_pressure"})
    gmgSmoother:set_forward_iter(gmgMainInverse)
    gmgSmoother:set_schur_iter(gmgPressureInverse)
    gmgSmoother:set_backward_iter(gmgMainInverse)

    local gmgSchurMassEq = ConvectionDiffusionFV1("p_pressure", "Lumen")
    gmgSchurMassEq:set_mass_scale(1.0)
    gmgSchurMassEq:set_diffusion(0.0)
    local gmgSchurMassDisc = DomainDiscretization(space)
    gmgSchurMassDisc:add(gmgSchurMassEq)
    local gmgSchurMembrane = ConvectionDiffusionFV1(
        "p_pressure", "Membrane"
    )
    gmgSchurMembrane:set_mass_scale(0.0)
    gmgSchurMembrane:set_diffusion(membraneDarcyMobility)
    local gmgSchurInter = ConvectionDiffusionFV1("p_pressure", "Inter")
    gmgSchurInter:set_mass_scale(0.0)
    gmgSchurInter:set_diffusion(interDarcyMobility)
    gmgSchurMassDisc:add(gmgSchurMembrane)
    gmgSchurMassDisc:add(gmgSchurInter)
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
            reduction = 1.0e-4,
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
    local componentNames = {"u", "v", "w", "p_pressure"}
    local function trueComponentDefects()
        solution:set_consistent_storage_type()
        local defect = rhs:clone()
        VecAssign(defect, rhs)
        matrix:apply_sub(defect, solution)
        local norms = {}
        for selected = 1, #componentNames do
            local factors = {0.0, 0.0, 0.0, 0.0}
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

-- The shared unknown already is physical kinematic pressure.
physicalState = coupledState:clone()
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
    physicalState, "p_pressure", "Inlet"
) / Integral(1.0, physicalState, "Inlet")
local outletPressureMean = Integral(
    physicalState, "p_pressure", "Outlet"
) / Integral(1.0, physicalState, "Outlet")
print("Apical pressure-trace jump: exactly zero (one shared pressure DoF)")
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
stokesOut:select_nodal("p_pressure", "stokes_pressure_kinematic")
stokesOut:print_subsets(
    vtk_file_name .. "_lumen", physicalState, "Lumen", 0, 0.0
)

local function writeDarcyOutput(suffix, subset, velocity, components)
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal("p_pressure", "darcy_pressure_kinematic")
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
