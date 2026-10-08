-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Lumen-flow test for the projector-backed 3rd_trace_new_build_test mesh.
--
-- The coarse UGX is loaded deliberately.  -numRefs constructs the conforming
-- projected hierarchy in memory, which is required by geometric multigrid.
-- Default units are consistent with the generated mesh: mm and seconds.

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
AssertPluginsLoaded({"NavierStokes", "SuperLU6", "neuro_collection"})

local totalStart = GetClockS()
local dim = 3
-- The whole-domain mesh remains runs/3rd_trace_new_build_test.ugx.  Flow is
-- solved on its Lua-extracted, regular lumen companion because the membrane's
-- constrained interface objects are outside the fluid domain and are not
-- accepted by NavierStokesFV1's regular-grid assembly.
local gridName = util.GetParam("-grid", "runs/3rd_trace_new_build_test_lumen.ugx")
local numRefs = math.floor(util.GetParamNumber("-numRefs", 1))
local numPreRefs = 0
local viscosity = util.GetParamNumber("-visc", 1.0)       -- mm^2/s
local inflow = util.GetParamNumber("-inflow", 8.4)        -- mm/s
local inletRadius = util.GetParamNumber("-inletRadius", 0.001) -- mm
local writeVTK = util.GetParamNumber("-writeVTK", 1) ~= 0
local upwind = util.GetParam("-upwind", "lps")
local stab = util.GetParam("-stab", "fields")
local diffLength = util.GetParam("-difflength", "cor")

assert(numRefs >= 1 and numRefs <= 2, "numRefs must be 1 or 2 for this test")
assert(viscosity > 0 and inletRadius > 0, "viscosity and inletRadius must be positive")

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print("Grid: " .. gridName)
print("In-memory projected refinements: " .. numRefs)
print("Equation: stationary Stokes")

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

local ns = NavierStokesFV1({"u", "v", "w", "p"}, {"Lumen"})
ns:set_exact_jacobian(false)
ns:set_stokes(true)
ns:set_laplace(true)
ns:set_kinematic_viscosity(viscosity)
ns:set_upwind(upwind)
ns:set_peclet_blend(false)
ns:set_stabilization(stab, diffLength)
ns:set_pac_upwind(false)

-- New-build inlet: center of the z-min terminal patch.  The terminal extender
-- creates a straight wall-normal tail, so the inward flow direction is +z.
local inletCenterX = 0.3248241776875
local inletCenterY = 0.46847242884375
local inletCenterZ = 0.2108

function newBuildInletVelocity(x, y, z, t)
    local dx, dy = x - inletCenterX, y - inletCenterY
    local r2 = dx * dx + dy * dy
    local q = r2 / (inletRadius * inletRadius)
    local profile = math.max(0.0, 1.0 - q)
    if q >= 1.0 - 1.0e-6 then profile = 0.0 end
    return 0.0, 0.0, inflow * profile
end

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

util.solver.defaults.approxSpace = approxSpace
local levelSmoother = ILU()
levelSmoother:set_damp(0.7)
levelSmoother:set_inversion_eps(1.0e-14)
levelSmoother:enable_overlap(true)
levelSmoother:enable_consistent_interfaces(true)

local gmg = GeometricMultiGrid(approxSpace)
gmg:set_base_level(numPreRefs)
gmg:set_base_solver(AgglomeratingSolver(SuperLU()))
gmg:set_gathered_base_solver_if_ambiguous(true)
gmg:set_smoother(levelSmoother)
gmg:set_cycle_type("V")
gmg:set_num_presmooth(4)
gmg:set_num_postsmooth(4)
gmg:set_rap(false)
gmg:set_transfer(StdTransfer())

local solver = util.solver.CreateSolver({
    type = "newton",
    convCheck = {
        type = "standard", iterations = 8,
        absolute = 1e-10, reduction = 1e-8, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = gmg,
        convCheck = {
            type = "standard", iterations = 200,
            absolute = 1e-10, reduction = 1e-8, verbose = true
        }
    }
})

local u = GridFunction(approxSpace)
for _, fct in ipairs({"u", "v", "w", "p"}) do Interpolate(0.0, u, fct) end

local gridStem = gridName:match("([^/\\]+)%.ugx$") or "3rd_trace_new_build_test"
local vtkPrefix = util.GetParam(
    "-vtk",
    MODEL_ROOT .. "/Results/new_build_" ..
        gridStem .. "_flow_ref" .. numRefs
)
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
