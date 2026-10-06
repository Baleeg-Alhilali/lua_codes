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

-- Standalone monolithic Stokes--Darcy--BJS model with an all-iterative
-- BiCGStab/GMG solver. No other application script is loaded by this file.
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
gridS = "long_tube"
gridName = "runs/" .. gridS .. ".ugx"
numPreRefs = 0
inletRadius = 1.0                 -- [mm]
tubeLength = 20.0                 -- [mm]
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
ilutThreshold = 1.0e-8
gmgBaseLevel = 0
gmgCycle = "V"
gmgPreSmooth = 1
gmgPostSmooth = 1
gmgUseRAP = true

-- MPI runs create .pvtu master files and one .vtu file per rank.
local flowModelTag = bstokes and "stokes" or "navier_stokes"
vtk_file_name = "Results/" .. gridS .. "_" .. flowModelTag
    .. "_darcy_bjs_gmg_ref" .. numRefs

assert(numRefs >= 0, "numRefs must be non-negative")
assert(vmax >= 0.0, "vmax must be non-negative")
assert(permeabilityM2 > 0.0, "permeabilityM2 must be positive")
assert(bjsAlpha > 0.0, "bjsAlpha must be positive")
assert(density > 0.0, "density must be positive")
assert(linearTolerance > 0.0 and linearReduction > 0.0,
    "linear convergence tolerances must be positive")
assert(linearMaxIters > 0 and coarseMaxIters > 0,
    "solver iteration limits must be positive")

-- Build every projected level before distributing the finished hierarchy.
balancer.firstDistLvl = numRefs

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
assert(file_contains(gridName, '<projector type="cylinder"'),
    "The long-tube grid has no serialized CylinderProjector: " .. gridName)

print("Cylinder-projector grid: " .. gridName)

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom,
    {"Lumen", "Inter", "Apical", "OuterWall", "Inlet", "Outlet",
     "TripleJunction"}),
    "The domain is missing a required Stokes--Darcy subset")

-- Repair open-axis vertices mislabeled Inter although they border Lumen.
local lumenSubset = dom:subset_handler():get_subset_index("Lumen")
local axisTolerance = 1.0e-10
local capTolerance = 1.0e-10
local function repair_lumen_axis_vertices()
    AssignSubset_VerticesInCube(
        dom,
        Vec3d(-axisTolerance, -axisTolerance, capTolerance),
        Vec3d( axisTolerance,  axisTolerance, tubeLength - capTolerance),
        lumenSubset
    )
end

repair_lumen_axis_vertices()
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)
repair_lumen_axis_vertices()

print("Domain info:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- Monolithic approximation space
--------------------------------------------------------------------------------

approxSpace = ApproximationSpace(dom)
local allSubsets = {
    "Lumen", "Inter", "Inlet", "Outlet", "Apical", "OuterWall",
    "TripleJunction"
}
local allSubsetsString =
    "Lumen,Inter,Inlet,Outlet,Apical,OuterWall,TripleJunction"

-- Fields exist on both adjacent elements so coupled exports can assemble
-- cross derivatives. Auxiliary extensions are pinned below.
approxSpace:add_fct(
    {"u", "v", "w", "pressure"}, "Lagrange", 1, allSubsets
)
approxSpace:add_fct(
    "darcyPressure", "Lagrange", 1, allSubsetsString
)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Physical coefficients and boundary data
--------------------------------------------------------------------------------

local permeability = permeabilityM2 * 1.0e6 -- [mm^2]
local darcyMobility = permeability / viscosity
local bjsFriction = bjsAlpha * viscosity / math.sqrt(permeability)
local outerPressure = outerPressurePa / (density * 1.0e-6)

print(string.format(
    "K=%.6e mm^2, K/nu=%.6e s, BJS friction=%.6e mm/s",
    permeability, darcyMobility, bjsFriction
))

function inletVelocity(x, y, z, t)
    local radiusSquared = x * x + y * y
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return 0.0, 0.0, vmax * profile
end

local function radialNormal(x, y)
    local radius = math.sqrt(x * x + y * y)
    if radius < 1.0e-12 then return 1.0, 0.0 end
    return x / radius, y / radius
end

function lumenNormalX(x, y, z, t, si)
    local nx = radialNormal(x, y)
    return nx
end
function lumenNormalY(x, y, z, t, si)
    local _, ny = radialNormal(x, y)
    return ny
end
function lumenNormalZ(x, y, z, t, si)
    return 0.0
end
function darcyNormalX(x, y, z, t, si)
    return -lumenNormalX(x, y, z, t, si)
end
function darcyNormalY(x, y, z, t, si)
    return -lumenNormalY(x, y, z, t, si)
end
function darcyNormalZ(x, y, z, t, si)
    return 0.0
end

--------------------------------------------------------------------------------
-- Stokes/Navier--Stokes problem in Lumen
--------------------------------------------------------------------------------

stokesEq = NavierStokesFV1({"u", "v", "w", "pressure"}, {"Lumen"})
-- Always compute the Stokes solution first.  For bstokes=false it is the
-- initial guess for the nonlinear Navier--Stokes correction below.
stokesEq:set_stokes(true)
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
stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

-- BJS--Saffman: (sigma_S n)_t = -alpha*nu/sqrt(K) u_S,t.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)

--------------------------------------------------------------------------------
-- Darcy problem and monolithic interface coupling
--------------------------------------------------------------------------------

darcyEq = ConvectionDiffusionFV1("darcyPressure", "Inter")
darcyEq:set_mass_scale(0.0)
darcyEq:set_diffusion(darcyMobility)

coupledDisc = DomainDiscretization(approxSpace)
coupledDisc:add(stokesEq)
coupledDisc:add(stokesInlet)
coupledDisc:add(stokesOpenBoundary)
coupledDisc:add(stokesRim)
coupledDisc:add(stokesBJS)

-- Normal-stress balance: sigma_S n_S = -p_D n_S.
local darcyPressureValue = darcyEq:value()
for _, row in ipairs({
    {"u", LuaUserNumber3d(lumenNormalX)},
    {"v", LuaUserNumber3d(lumenNormalY)},
    {"w", LuaUserNumber3d(lumenNormalZ)}
}) do
    local interfaceTraction = NeumannBoundary(row[1], "fv1")
    interfaceTraction:add(darcyPressureValue * row[2], "Apical", "Lumen")
    coupledDisc:add(interfaceTraction)
end

-- Normal-flux continuity using opposite outward interface orientations.
local stokesVelocity = stokesEq:velocity()
local velocityComponents = {}
for component = 0, 2 do
    velocityComponents[component + 1] = UserVectorEntryAdapter3d()
    velocityComponents[component + 1]:set_vector(stokesVelocity, component)
end
local darcyNormalVelocity =
      velocityComponents[1] * LuaUserNumber3d(darcyNormalX)
    + velocityComponents[2] * LuaUserNumber3d(darcyNormalY)
    + velocityComponents[3] * LuaUserNumber3d(darcyNormalZ)

darcyInterfaceFlux = NeumannBoundary("darcyPressure", "fv1")
darcyInterfaceFlux:add(darcyNormalVelocity, "Apical", "Inter")
darcyOuterPressure = DirichletBoundary()
darcyOuterPressure:add(outerPressure, "darcyPressure", "OuterWall")

inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(0.0, "darcyPressure", "Lumen,Inlet,Outlet")
inactiveStokes = DirichletBoundary()
for _, fct in ipairs({"u", "v", "w", "pressure"}) do
    inactiveStokes:add(0.0, fct, "Inter,OuterWall")
end

coupledDisc:add(darcyEq)
coupledDisc:add(darcyInterfaceFlux)
coupledDisc:add(darcyOuterPressure)
coupledDisc:add(inactiveDarcy)
coupledDisc:add(inactiveStokes)

--------------------------------------------------------------------------------
-- BiCGStab with geometric multigrid preconditioning
--------------------------------------------------------------------------------

local solve_start_time = GetClockS()
util.solver.defaults.approxSpace = approxSpace
local coupledMatrix = AssembledLinearOperator(coupledDisc)
local coupledState = GridFunction(approxSpace)
local coupledRhs = GridFunction(approxSpace)
coupledState:set(0.0)
coupledDisc:adjust_solution(coupledState)
coupledDisc:assemble_linear(coupledMatrix, coupledRhs)

local gmgPreconditioner = {
    type = "gmg",
    approxSpace = approxSpace,
    discretization = coupledDisc,
    baseLevel = gmgBaseLevel,
    baseSolver = {
        type = "bicgstab",
        precond = {
            type = "ilut",
            threshold = ilutThreshold,
            ordering = "NativeCuthillMcKee"
        },
        convCheck = {
            type = "standard",
            iterations = coarseMaxIters,
            absolute = coarseTolerance,
            reduction = coarseReduction,
            verbose = true
        }
    },
    smoother = {
        type = "ilut",
        threshold = ilutThreshold,
        ordering = "NativeCuthillMcKee"
    },
    cycle = gmgCycle,
    preSmooth = gmgPreSmooth,
    postSmooth = gmgPostSmooth,
    rap = gmgUseRAP,
    transfer = {
        type = "std",
        restrictionDamp = 1.0,
        prolongationDamp = 1.0,
        enableP1LagrangeOptimization = true
    }
}

local function outerLinearSolverDescription()
    return {
        type = "bicgstab",
        precond = gmgPreconditioner,
        convCheck = {
            type = "standard",
            iterations = linearMaxIters,
            absolute = linearTolerance,
            reduction = linearReduction,
            verbose = true
        }
    }
end

print("Solving the monolithic Stokes--Darcy initial field ...")
local iterativeSolver = util.solver.CreateLinearSolver(
    outerLinearSolverDescription()
)

iterativeSolver:init(coupledMatrix, coupledState)
assert(iterativeSolver:apply(coupledState, coupledRhs),
    "BiCGStab--GMG Stokes--Darcy solver did not converge")

if not bstokes then
    stokesEq:set_stokes(false)
    local nonlinearSolver = util.solver.CreateSolver({
        type = "newton",
        linSolver = outerLinearSolverDescription(),
        lineSearch = {
            type = "standard",
            maxSteps = 10,
            lambdaStart = 1.0,
            lambdaReduce = 0.5,
            acceptBest = true,
            checkAll = false
        },
        convCheck = {
            type = "standard",
            iterations = newtonMaxIters,
            absolute = newtonTolerance,
            reduction = newtonReduction,
            verbose = true
        }
    })
    nonlinearSolver:init(AssembledOperator(coupledDisc))
    nonlinearSolver:prepare(coupledState)
    print("Applying the monolithic Navier--Stokes--Darcy correction ...")
    assert(nonlinearSolver:apply(coupledState),
        "Newton/BiCGStab--GMG Navier--Stokes--Darcy solver did not converge")
end

--------------------------------------------------------------------------------
-- Darcy velocity and VTK output
--------------------------------------------------------------------------------

darcyVelocity = DarcyVelocityLinker()
darcyVelocity:set_permeability(permeability)
darcyVelocity:set_viscosity(viscosity)
darcyVelocity:set_density(1.0)
darcyVelocity:set_gravity(ConstUserVector(0.0))
darcyVelocity:set_pressure_gradient(darcyEq:gradient())

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
    vtk_file_name .. "_lumen", coupledState, "Lumen", 0, 0.0
)

darcyOut = VTKOutput()
darcyOut:clear_selection()
darcyOut:select_nodal("darcyPressure", "darcy_pressure_kinematic")
darcyOut:select_element(darcyVelocity, "darcy_velocity_mm_s")
darcyOut:select_element(darcyVelocityComponents[1], "u")
darcyOut:select_element(darcyVelocityComponents[2], "v")
darcyOut:select_element(darcyVelocityComponents[3], "w")
darcyOut:print_subsets(
    vtk_file_name .. "_inter", coupledState, "Inter", 0, 0.0
)

local solve_elapsed_time = GetClockS() - solve_start_time
local total_elapsed_time = GetClockS() - total_start_time
print(string.format(
    "Coupled Stokes--Darcy GMG runtime: %.6f mins",
    solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:             %.6f mins",
    total_elapsed_time / 60.0
))
