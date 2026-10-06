-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- Standalone curved-tube Robin-reduced Stokes--Darcy--BJS model with
-- all-iterative BiCGStab/GMG subsolvers. No other application script is loaded.
AssertPluginsLoaded({"NavierStokes", "ConvectionDiffusion"})

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- User configuration
--------------------------------------------------------------------------------

dim = 3

-- These are the only command-line parameters accepted by this script.
numRefs = util.GetParamNumber("-numRefs", 2, "number of uniform refinements")
bstokes = util.GetParamBool("-bstokes", true,
    "true: Stokes flow; false: Navier--Stokes flow")
vmax = util.GetParamNumber("-vmax", 3.0, "maximum inlet velocity [mm/s]")

-- Geometry and physical parameters: edit values here, not on the command line.
gridS = "curved_tube"
gridName = "runs/" .. gridS .. ".ugx"
numPreRefs = 0
inletRadius = 1.0e-3              -- [mm], lumen radius measured on Apical
inletCenterY = -6.774645191483331e-4 -- [mm]
inletCenterZ = 1.949999995936047e-3  -- [mm]
tubeLength = 8.3281559e-2         -- approximate centerline length [mm]
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
permeabilityM2 = 1.0e-12          -- Darcy permeability [m^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- [kg/m^3]
outerPressurePa = 0.0             -- pressure at OuterWall [Pa]

-- Solver parameters: edit values here, not on the command line.
linearTolerance = 1.0e-10
linearReduction = 1.0e-8
linearMaxIters = 1000
newtonTolerance = 1.0e-10
newtonReduction = 1.0e-8
newtonMaxIters = 30
coarseTolerance = 1.0e-12
coarseReduction = 1.0e-10
coarseMaxIters = 500
ilutThreshold = 1.0e-12
coarseILUTInversionEps = 1.0e-20
gmgBaseLevel = 0
gmgCycle = "V"
gmgPreSmooth = 4
gmgPostSmooth = 4
gmgUseRAP = false

-- Porous-layer impedance parameter.
porousResistanceLength = 3.0e-2 -- [mm], approximate porous-layer thickness

-- MPI runs create .pvtu master files and one .vtu file per rank.
local flowModelTag = bstokes and "stokes" or "navier_stokes"
vtk_file_name = "Results/" .. gridS .. "_" .. flowModelTag
    .. "_darcy_bjs_robin_gmg_ref" .. numRefs

assert(numRefs >= 0, "numRefs must be non-negative")
assert(vmax >= 0.0, "vmax must be non-negative")
assert(permeabilityM2 > 0.0, "permeabilityM2 must be positive")
assert(bjsAlpha > 0.0, "bjsAlpha must be positive")
assert(density > 0.0, "density must be positive")
assert(linearTolerance > 0.0 and linearReduction > 0.0,
    "linear convergence tolerances must be positive")
assert(linearMaxIters > 0 and coarseMaxIters > 0,
    "solver iteration limits must be positive")
assert(porousResistanceLength > 0.0,
    "porousResistanceLength must be positive")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print(string.format(
    "Model: %s, vmax=%.6e mm/s, refinements=%d",
    bstokes and "Stokes" or "Navier--Stokes", vmax, numRefs
))

--------------------------------------------------------------------------------
-- Domain
--------------------------------------------------------------------------------

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
    "The curved grid has no serialized NeuriteProjector: " .. gridName)
assert(file_contains(gridName, '<vertex_attachment name="npSurfParams"'),
    "The curved grid has no npSurfParams refinement data: " .. gridName)

print("Neurite-projector curved grid: " .. gridName)

-- The first load registers the custom npSurfParams attachment type. The
-- second load can then restore its values and the serialized projector.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Membrane", "Inter", "Apical", "Basolateral",
     "OuterWall", "Inlet", "Outlet"}),
    "The domain is missing a required Stokes--Darcy subset")

-- Keep level 0 on one rank until the complete projected hierarchy exists.
-- Distributing the coarse projector grid first produced inconsistent MPI
-- refinement in the validated new-build flow test.
if NumProcs() > 1 then
    balancer.firstDistLvl = numRefs
    balancer.ParseParameters()
    print("First distribution level: " .. balancer.firstDistLvl)
end
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- Separate approximation spaces for a stable partitioned interface solve
--------------------------------------------------------------------------------

stokesSpace = ApproximationSpace(dom)
stokesSpace:add_fct(
    {"u", "v", "w", "pressure"}, "Lagrange", 1,
    {"Lumen", "Inlet", "Outlet", "Apical", "TripleJunction"}
)
stokesSpace:init_levels()
stokesSpace:init_top_surface()

darcySpace = ApproximationSpace(dom)
darcySpace:add_fct(
    "darcyPressure", "Lagrange", 1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet,TripleJunction"
)
darcySpace:init_levels()
darcySpace:init_top_surface()

print("Stokes approximation space:")
stokesSpace:print_statistic()
print("Darcy approximation space:")
darcySpace:print_statistic()

--------------------------------------------------------------------------------
-- Physical coefficients and boundary data
--------------------------------------------------------------------------------

local permeability = permeabilityM2 * 1.0e6 -- [mm^2]
local darcyMobility = permeability / viscosity
local bjsFriction = bjsAlpha * viscosity / math.sqrt(permeability)
local outerPressure = outerPressurePa / (density * 1.0e-6)
local normalRobinResistance = porousResistanceLength / darcyMobility

print(string.format(
    "K=%.6e mm^2, K/nu=%.6e s, BJS friction=%.6e mm/s, Robin=%.6e mm/s",
    permeability, darcyMobility, bjsFriction, normalRobinResistance
))

function inletVelocity(x, y, z, t)
    -- The inlet is the x-min cap. Its inward tangent is +x.
    local dy, dz = y - inletCenterY, z - inletCenterZ
    local radiusSquared = dy * dy + dz * dz
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return vmax * profile, 0.0, 0.0
end

--------------------------------------------------------------------------------
-- One-way interface data: Stokes normal velocity drives Darcy flux.
--------------------------------------------------------------------------------

stokesState = GridFunction(stokesSpace)
darcyState = GridFunction(darcySpace)
stokesState:set(0.0)
darcyState:set(0.0)

local transferredDarcyFlux = 0.0

function laggedStokesNormalFlux(x, y, z, t, si)
    return transferredDarcyFlux
end

--------------------------------------------------------------------------------
-- Stokes/Navier--Stokes problem in Lumen
--------------------------------------------------------------------------------

stokesEq = NavierStokesFV1({"u", "v", "w", "pressure"}, {"Lumen"})
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
stokesOpenBoundary:add("Outlet,Apical")

-- BJS--Saffman: (sigma_S n)_t = -alpha*nu/sqrt(K) u_S,t.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalRobinResistance)

-- Close the cap/wall junction explicitly. TripleJunction is an edge subset in
-- this 3D mesh; NavierStokesWall constrains its velocity DoFs so they are not
-- left to the competing inlet/outlet and Apical interface conditions.
stokesTripleJunctionWall = NavierStokesWall(stokesEq)
stokesTripleJunctionWall:add("TripleJunction")

-- Reduced porous-layer condition in the interface-normal direction:
--     (sigma_S n).n + R_D u_n = 0,
-- where R_D = L_D/(K/nu).  The WSBC operator assembles this term into both
-- defect and Jacobian, unlike a prescribed Neumann callback.

--------------------------------------------------------------------------------
-- Apply the normal Darcy impedance implicitly; BJS handles tangential slip.
stokesDisc = DomainDiscretization(stokesSpace)
stokesDisc:add(stokesEq)
stokesDisc:add(stokesInlet)
stokesDisc:add(stokesOpenBoundary)
stokesDisc:add(stokesBJS)
stokesDisc:add(stokesTripleJunctionWall)

--------------------------------------------------------------------------------
-- Darcy reconstruction driven by the stabilized Stokes normal flux.
--------------------------------------------------------------------------------

darcyMembraneEq = ConvectionDiffusionFV1("darcyPressure", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(darcyMobility)

darcyInterEq = ConvectionDiffusionFV1("darcyPressure", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(darcyMobility)

darcyInterfaceFlux = NeumannBoundary("darcyPressure", "fv1")
darcyInterfaceFlux:add("laggedStokesNormalFlux", "Apical", "Membrane")
darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(outerPressure, "darcyPressure", "OuterWall")

darcyDisc = DomainDiscretization(darcySpace)
darcyDisc:add(darcyMembraneEq)
darcyDisc:add(darcyInterEq)
darcyDisc:add(darcyInterfaceFlux)
darcyDisc:add(darcyOuterPressure)
--------------------------------------------------------------------------------
-- Sequential BiCGStab/GMG solves
--------------------------------------------------------------------------------

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    -- The level-0 matrix is distributed as well.  Merely asking GMG to
    -- "gather if ambiguous" left BiCGStab working on an incomplete parallel
    -- coarse problem and it stagnated at numRefs=2.  Build the coarse inverse
    -- explicitly on the agglomerated matrix.  Its actual solve remains fully
    -- iterative (BiCGStab/ILUT); no LU factorization is used.
    local iterativeCoarseSolver = util.solver.CreateLinearSolver({
        type = "bicgstab",
        precond = {
            type = "ilut",
            threshold = ilutThreshold,
            ordering = "NativeCuthillMcKee",
            inversionEps = coarseILUTInversionEps
        },
        convCheck = {
            type = "standard",
            iterations = coarseMaxIters,
            absolute = coarseTolerance,
            reduction = coarseReduction,
            verbose = true
        }
    })
    local agglomeratedIterativeCoarseSolver =
        AgglomeratingSolver(iterativeCoarseSolver)

    return {
        type = "bicgstab",
        precond = {
            type = "gmg",
            approxSpace = space,
            discretization = disc,
            baseLevel = gmgBaseLevel,
            baseSolver = agglomeratedIterativeCoarseSolver,
            smoother = {
                type = "ilu",
                damping = 0.7,
                overlap = true,
                consistentInterfaces = true,
                inversionEps = 1.0e-12
            },
            cycle = gmgCycle,
            preSmooth = gmgPreSmooth,
            postSmooth = gmgPostSmooth,
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
            iterations = linearMaxIters,
            absolute = linearTolerance,
            reduction = linearReduction,
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
    local solver = util.solver.CreateLinearSolver(
        outerLinearSolverDescription(space, disc)
    )
    print("=== " .. label .. " OUTER BiCGStab/GMG convergence ===")
    solver:init(matrix, solution)
    assert(solver:apply(solution, rhs), label .. " BiCGStab--GMG failed")
    print("=== " .. label .. " OUTER solve converged ===")
    solution:set_consistent_storage_type()
    return solution
end

for _, state in ipairs({stokesState, darcyState}) do
    state:set_consistent_storage_type()
end

local solvedStokes = solveLinearSubproblem(stokesSpace, stokesDisc, "Stokes")
VecScaleAdd2(stokesState, 1.0, solvedStokes, 0.0, stokesState)
stokesState:set_consistent_storage_type()

-- Conservative MPI-safe transfer. The integral is globally reduced over all
-- ranks, and its mean is imposed around the complete cylindrical interface.
local stokesVelocityData = GridFunctionVectorData(stokesState, "u,v,w")
local totalStokesOutflow = IntegralNormalComponentOnManifold(
    stokesVelocityData, stokesState, "Apical", "Lumen", 0.0
)
local apicalArea = 2.0 * math.pi * inletRadius * tubeLength
transferredDarcyFlux = -totalStokesOutflow / apicalArea
print(string.format(
    "Interface flux transfer: total=%.6e mm^3/s, mean Darcy-normal=%.6e mm/s",
    totalStokesOutflow, transferredDarcyFlux
))

local solvedDarcy = solveLinearSubproblem(darcySpace, darcyDisc, "Darcy")
VecScaleAdd2(darcyState, 1.0, solvedDarcy, 0.0, darcyState)
darcyState:set_consistent_storage_type()

--------------------------------------------------------------------------------
-- Darcy velocity and VTK output
--------------------------------------------------------------------------------

darcyVelocity = DarcyVelocityLinker()
darcyVelocity:set_permeability(permeability)
darcyVelocity:set_viscosity(viscosity)
darcyVelocity:set_density(1.0)
darcyVelocity:set_gravity(ConstUserVector(0.0))
-- Both volume discretizations use the same pressure field and coefficient.
-- The gradient data is identical in form on either side of Basolateral.
darcyVelocity:set_pressure_gradient(darcyInterEq:gradient())

local darcyVelocityComponents = {}
for component = 0, 2 do
    darcyVelocityComponents[component + 1] = UserVectorEntryAdapter3d()
    darcyVelocityComponents[component + 1]:set_vector(
        darcyVelocity, component
    )
end

stokesOut = VTKOutput()
stokesOut:clear_selection()
stokesOut:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
stokesOut:select_nodal("u", "u")
stokesOut:select_nodal("v", "v")
stokesOut:select_nodal("w", "w")
stokesOut:select_nodal("pressure", "stokes_pressure_kinematic")
stokesOut:print_subsets(
    vtk_file_name .. "_lumen", stokesState, "Lumen", 0, 0.0
)

darcyOut = VTKOutput()
darcyOut:clear_selection()
darcyOut:select_nodal("darcyPressure", "darcy_pressure_kinematic")
darcyOut:select_element(darcyVelocity, "darcy_velocity_mm_s")
darcyOut:select_element(darcyVelocityComponents[1], "u")
darcyOut:select_element(darcyVelocityComponents[2], "v")
darcyOut:select_element(darcyVelocityComponents[3], "w")
darcyOut:print_subsets(
    vtk_file_name .. "_inter", darcyState, "Membrane,Inter", 0, 0.0
)

local solve_elapsed_time = GetClockS() - solve_start_time
local total_elapsed_time = GetClockS() - total_start_time
print(string.format(
    "Partitioned Stokes--Darcy GMG runtime: %.6f mins",
    solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:             %.6f mins",
    total_elapsed_time / 60.0
))
