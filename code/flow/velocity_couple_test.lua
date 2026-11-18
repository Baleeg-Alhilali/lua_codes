--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Legacy monolithic Stokes--Darcy--BJS test on the cylinder-projector long
-- tube. It couples lumen traction and normal velocity to interstitial Darcy
-- pressure and uses an ILUT-preconditioned iterative linear solve.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolves model-relative paths and loads UG4 flow, diffusion, and balancing.
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

-- The coupled Stokes--Darcy operator is nonsymmetric and contains the Stokes
-- saddle-point block.  Solve it with preconditioned BiCGStab; no direct-solver
-- plugin is required.
AssertPluginsLoaded({"NavierStokes", "ConvectionDiffusion"})

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------
-- Only -numRefs and -output are command-line inputs. Edit the remaining
-- geometry, physics, discretization, and solver values in this block.

dim = 3

-- Always refine the original UGX written by the mesh-generation script.  This
-- straight test grid carries a serialized CylinderProjector (not a
-- NeuriteProjector).
gridS ="long_tube" 
gridName = "runs/"..gridS..".ugx"

numRefs = util.GetParamNumber("-numRefs", 2)
numPreRefs = 0 -- 0

-- Build every projected level before distributing the finished hierarchy.
-- This is the refinement order validated by the projector-flow tests.
balancer.firstDistLvl = numRefs

v_avg = 1.5 --mm/s
v_max = 2*v_avg --mm/s
-- Physical parameters
viscosity 	= 1.0               -- kinematic viscosity [mm^2/s]
inflow		= v_max             -- maximum inlet velocity [mm/s]
inletRadius = 1.0                -- inlet lumen radius [mm]
permeabilityM2 = 1.0e-12         -- Darcy permeability [m^2]
permeability = permeabilityM2 * 1.0e6 -- [mm^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- fluid density [kg/m^3]
outerPressurePa = 0.0             -- OuterWall pressure [Pa]
assert(permeabilityM2 > 0.0, "permeabilityM2 must be positive")
assert(bjsAlpha > 0.0, "bjsAlpha must be positive")
assert(density > 0.0, "density must be positive")
vtk_file_name = util.GetParam("-output",
    "Results/" .. gridS .. "_stokes_darcy_bjs_ref" .. numRefs,
    "VTK output prefix")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Loads/refines the cylinder-projector mesh, repairs its legacy lumen-axis
-- labels, and defines one space for Stokes and Darcy unknowns.

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
assert(util.CheckSubsets(dom, {"Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"}),
       "The domain is missing a required lumen subset")

-- long_tube.ugx labels its interior axis vertices as Inter even though they
-- are vertices of the central Lumen prisms.  A subset-restricted
-- approximation space then has no flow DoF at one corner of every such
-- prism.  Repair only the open axis; keep the two cap-centre vertices in
-- Inlet and Outlet.
local lumenSubset = dom:subset_handler():get_subset_index("Lumen")
local axisTolerance = 1.0e-10
local capTolerance = 1.0e-10
local tubeLength = 20.0
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

approxSpace = ApproximationSpace(dom)
local allSubsets = {
    "Lumen", "Inter", "Inlet", "Outlet", "Apical", "OuterWall",
    "TripleJunction"
}
local allSubsetsString =
    "Lumen,Inter,Inlet,Outlet,Apical,OuterWall,TripleJunction"
approxSpace:add_fct(
    {"u", "v", "w", "pressure"}, "Lagrange", 1,
    allSubsets
)
approxSpace:add_fct(
    "darcyPressure", "Lagrange", 1,
    allSubsetsString
)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

-- Order the DoFs:
-- OrderLex crashes in DoFDistribution::permute_indices for this approximation
-- space in the selected UG4 build.  The default ordering is valid.
-- OrderLex(approxSpace, "xy")

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Defines all inlet/normal functions first, converts physical coefficients,
-- and then assembles Stokes, Darcy, BJS, traction, and flux-continuity terms.

-- Pressures are kinematic (mm^2/s^2).  Darcy's law is
-- u_D = -(K/nu) grad(p_D).  The BJS--Saffman coefficient below has units mm/s.
local darcyMobility = permeability / viscosity
local bjsFriction = bjsAlpha * viscosity / math.sqrt(permeability)
local outerPressure = outerPressurePa / (density * 1.0e-6)

print(string.format(
    "K=%.6e mm^2, K/nu=%.6e s, BJS friction=%.6e mm/s",
    permeability, darcyMobility, bjsFriction
))

function inletVelocity(x, y, z, t)
    local radiusSquared = x * x + y * y
    local profile = math.max(0.0, 1.0 - radiusSquared / (inletRadius * inletRadius))
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return 0.0, 0.0, inflow * profile
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
-- Stokes equation and interface conditions
--------------------------------------------------------------------------------

stokesEq = NavierStokesFV1({"u", "v", "w", "pressure"}, {"Lumen"})
stokesEq:set_stokes(true)
stokesEq:set_laplace(false) -- symmetric-gradient stress is required for BJS
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_stabilization("fields", "cor")

stokesInlet = NavierStokesInflow(stokesEq)
stokesInlet:add("inletVelocity", "Inlet")
stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
stokesOpenBoundary:add("Outlet,Apical")
stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")
-- BJS--Saffman: (sigma_S n)_t = -alpha*nu/sqrt(K) * u_S,t.
-- NavierStokesWSBCFV1 assembles the corresponding positive tangential
-- dissipation in the Stokes operator, so this stiff term is treated
-- implicitly instead of being lagged in the interface iteration.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(bjsFriction)
stokesBJS:set_sliding_limit(0.0)

darcyEq = ConvectionDiffusionFV1("darcyPressure", "Inter")
darcyEq:set_mass_scale(0.0)
darcyEq:set_diffusion(darcyMobility)

coupledDisc = DomainDiscretization(approxSpace)
coupledDisc:add(stokesEq)
coupledDisc:add(stokesInlet)
coupledDisc:add(stokesOpenBoundary)
coupledDisc:add(stokesRim)
coupledDisc:add(stokesBJS)

-- Normal-stress balance: sigma_S n_S = -p_D n_S.  NeumannBoundary
-- subtracts the supplied outward PDE flux, so p_D e_i gives the required
-- component p_D n_i while retaining the derivative with respect to p_D.
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

--------------------------------------------------------------------------------
-- Darcy equation through Inter
--------------------------------------------------------------------------------

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

-- The coupled data exports require each field on both adjacent elements.
-- Pin the auxiliary extensions away from Apical; the physical equations still
-- act only on Lumen (Stokes) and Inter (Darcy).
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
-- Solver
--------------------------------------------------------------------------------
-- Assembles the monolithic linear system and solves it with BiCGStab/ILUT.

local solve_start_time = GetClockS()
util.solver.defaults.approxSpace = approxSpace
local coupledMatrix = AssembledLinearOperator(coupledDisc)
local coupledState = GridFunction(approxSpace)
local coupledRhs = GridFunction(approxSpace)
coupledState:set(0.0)
coupledDisc:adjust_solution(coupledState)
coupledDisc:assemble_linear(coupledMatrix, coupledRhs)

local solverDesc = {
    type = "bicgstab", -- alternatives: "gmres" or "cg" for compatible operators
    precond = {
        type = "ilut", -- alternative: "ilu"; use the dedicated GMG script for GMG
        threshold = 1.0e-8,
        ordering = "NativeCuthillMcKee"
    },
    convCheck = {
        type = "standard",
        iterations = 1000,
        absolute = 1e-14,
        reduction = 1.0e-8,
        verbose = true
    }
}
local iterativeSolver = util.solver.CreateLinearSolver(solverDesc)
iterativeSolver:init(coupledMatrix, coupledState)
assert(iterativeSolver:apply(coupledState, coupledRhs),
    "Iterative Stokes--Darcy solver did not converge")

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- This is a steady problem. No time discretization is used; coupledState is
-- initialized to zero before the linear system is assembled.

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- Reconstructs Darcy velocity and writes separate lumen and interstitial VTK.

darcyVelocity = DarcyVelocityLinker()
darcyVelocity:set_permeability(permeability)
darcyVelocity:set_viscosity(viscosity)
darcyVelocity:set_density(1.0)
darcyVelocity:set_gravity(ConstUserVector(0.0))
darcyVelocity:set_pressure_gradient(darcyEq:gradient())

stokesOut = VTKOutput()
stokesOut:clear_selection()
stokesOut:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
stokesOut:select_nodal("pressure", "stokes_pressure_kinematic")
stokesOut:print_subsets(vtk_file_name .. "_lumen", coupledState, "Lumen", 0, 0.0)

darcyOut = VTKOutput()
darcyOut:clear_selection()
darcyOut:select_nodal("darcyPressure", "darcy_pressure_kinematic")
darcyOut:select_element(darcyVelocity, "darcy_velocity_mm_s")
darcyOut:print_subsets(vtk_file_name .. "_inter", coupledState, "Inter", 0, 0.0)

local solve_elapsed_time = GetClockS() - solve_start_time
local total_elapsed_time = GetClockS() - total_start_time
print(string.format("Coupled Stokes--Darcy runtime: %.6f mins", solve_elapsed_time / 60.0))
print(string.format("Total script runtime:        %.6f mins", total_elapsed_time / 60.0))
