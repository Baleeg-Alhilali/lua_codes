--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Solves steady three-dimensional Stokes flow in a projector-backed nephron
-- lumen using FV1 stabilization and a geometric-multigrid linear solver.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolves model-relative paths and loads UG4 utilities and required plugins.
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
AssertPluginsLoaded({"NavierStokes", "SuperLU6", "neuro_collection"})

local totalStart = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------
-- Only -numRefs and -output are command-line inputs. All physics, mesh,
-- discretization, and solver choices are maintained here.
local dim = 3
-- The whole-domain mesh remains runs/3rd_trace_new_build_test.ugx.  Flow is
-- solved on its Lua-extracted, regular lumen companion because the membrane's
-- constrained interface objects are outside the fluid domain and are not
-- accepted by NavierStokesFV1's regular-grid assembly.
local gridName = "runs/third_trace_third_test2.ugx"
local numRefs = math.floor(util.GetParamNumber("-numRefs", 2))
local numPreRefs = 0
local viscosity = 1.0              -- [mm^2/s]
local inflow = 8.4                 -- maximum inlet velocity [mm/s]
local inletRadius = 0.015          -- [mm]
local writeVTK = true
local upwind = "lps"
local stab = "fields"
local diffLength = "cor"
local outputFileName = util.GetParam(
    "-output",
    MODEL_ROOT .. "/Results_final/new_build_flow_ref" .. numRefs,
    "VTK output prefix"
)

assert(numRefs >= 1 and numRefs <= 2, "numRefs must be 1 or 2 for this test")
assert(viscosity > 0 and inletRadius > 0, "viscosity and inletRadius must be positive")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print("Grid: " .. gridName)
print("In-memory projected refinements: " .. numRefs)
print("Equation: stationary Stokes")

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Loads the coarse projected grid, builds its in-memory refinement hierarchy,
-- and defines the lumen velocity-pressure approximation space.

local function fileContains(path, needle)
    local f = assert(io.open(path, "r"), "Cannot open grid file: " .. path)
    for line in f:lines() do
        if line:find(needle, 1, true) then f:close(); return true end
    end
    f:close()
    return false
end

assert(fileContains(gridName, '<projector type="neurite"'),
       "Input UGX has no serialized NeuriteProjector")
assert(fileContains(gridName, 'name="npSurfParams"'),
       "Input UGX has no npSurfParams attachment")

-- The first load registers the custom attachment type.  The second load then
-- restores all npSurfParams values and the serialized projector.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

local dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {"Lumen", "Apical", "Inlet", "Outlet"}),
       "Required lumen boundary subsets are missing")

-- Keep level 0 on one rank until projected level 1 exists.  Distribution of
-- the coarse grid before projection previously produced an inconsistent MPI
-- hierarchy.
if NumProcs() > 1 then
    balancer.firstDistLvl = numRefs
    balancer.ParseParameters()
    print("First distribution level: " .. balancer.firstDistLvl)
end
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)
print(dom:domain_info():to_string())

local approxSpace = ApproximationSpace(dom)
for _, fct in ipairs({"u", "v", "w", "p"}) do
    approxSpace:add_fct(fct, "Lagrange", 1, "Lumen,Inlet,Outlet,Apical,TripleJunction")
end
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Defines the inlet function first, then assembles Stokes, inflow, open-outlet,
-- and impermeable-wall discretizations.

-- Inlet of third_trace_third_test2.ugx. The cap lies in the z=1.127 plane,
-- and the adjacent lumen layers extend in +z, so inward flow points in +z.
local inletCenterX = 1.91647690475
local inletCenterY = 2.782189332375
local inletCenterZ = 1.127

function newBuildInletVelocity(x, y, z, t)
    local dx, dy = x - inletCenterX, y - inletCenterY
    local r2 = dx * dx + dy * dy
    local q = r2 / (inletRadius * inletRadius)
    local profile = math.max(0.0, 1.0 - q)
    if q >= 1.0 - 1.0e-6 then profile = 0.0 end
    return 0.0, 0.0, inflow * profile
end

local ns = NavierStokesFV1({"u", "v", "w", "p"}, {"Lumen"})
ns:set_exact_jacobian(false)
ns:set_stokes(true)
ns:set_laplace(true)
ns:set_kinematic_viscosity(viscosity)
ns:set_upwind(upwind)
ns:set_peclet_blend(false)
ns:set_stabilization(stab, diffLength)
ns:set_pac_upwind(false)

local inletDisc = NavierStokesInflow(ns)
inletDisc:add("newBuildInletVelocity", "Inlet")
local outletDisc = NavierStokesNoNormalStressOutflow(ns)
outletDisc:add("Outlet")
local wallDisc = NavierStokesWall(ns)
wallDisc:add("Apical,TripleJunction")

local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(ns)
domainDisc:add(inletDisc)
domainDisc:add(outletDisc)
domainDisc:add(wallDisc)

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- Configures Newton around the linear Stokes operator with BiCGStab/GMG and
-- a gathered SuperLU coarse solve.

util.solver.defaults.approxSpace = approxSpace
local solverDesc = {
    type = "newton", -- alternative for this linear Stokes case: a linear solver
    convCheck = {
        type = "standard",
        iterations = 8,
        absolute = 1e-14,
        reduction = 1.0e-8,
        verbose = true
    },
    linSolver = {
        type = "bicgstab", -- alternative: "gmres"
        precond = {
            type = "gmg", -- alternative: "ilu" for a single-level solve
            approxSpace = approxSpace,
            baseLevel = 0,
            baseSolver = "superlu", -- alternative: iterative coarse solver
            gatheredBaseSolverIfAmbiguous = true,
            smoother = {
                type = "ilu", -- alternative: "jac"
                damping = 0.7,
                inversionEps = 1.0e-14,
                overlap = true,
                consistentInterfaces = true
            },
            cycle = "V", -- alternative: "W"
            preSmooth = 4,
            postSmooth = 4,
            rap = false,
            transfer = "std"
        },
        convCheck = {
            type = "standard",
            iterations = 200,
            absolute = 1e-14,
            reduction = 1.0e-8,
            verbose = true
        }
    }
}
local solver = util.solver.CreateSolver(solverDesc)

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- This is a steady problem. No time discretization is used; every unknown is
-- initialized to zero before the stationary solve.

local u = GridFunction(approxSpace)
for _, fct in ipairs({"u", "v", "w", "p"}) do Interpolate(0.0, u, fct) end

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- Writes the initialized and converged lumen fields and reports solve timing.

local vtkPrefix = outputFileName
local out
if writeVTK then
    print("Flow output prefix: " .. vtkPrefix)
    out = VTKOutput()
    out:clear_selection()
    out:select_nodal({"u", "v", "w"}, "velocity")
    out:select_nodal("u", "u")
    out:select_nodal("v", "v")
    out:select_nodal("w", "w")
    out:select_nodal("p", "pressure")
    out:print_subsets(vtkPrefix .. "_lumen", u, "Lumen", 0, 0.0)
end

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)
local solveStart = GetClockS()
local ok = solver:apply(u)
assert(ok ~= false, "Stationary Stokes solve reported failure")

if writeVTK then
    out:print_subsets(vtkPrefix .. "_lumen", u, "Lumen", 1, 0.0)
    out:write_time_pvd(vtkPrefix .. "_lumen", u)
end

print(string.format("FLOW_TEST_SUCCESS refs=%d ranks=%d solve_s=%.6f total_s=%.6f",
    numRefs, NumProcs(), GetClockS() - solveStart, GetClockS() - totalStart))
