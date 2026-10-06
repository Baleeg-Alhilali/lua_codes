-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Working copy for Stokes--Darcy coupling development.
-- This script still solves only Lumen Stokes flow; Apical is a wall.
-- The partitioned_experimental.lua trial did not pass interface validation.
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- The level-0 Stokes/Navier--Stokes system is a saddle-point problem.  Its
-- scalar pressure rows are not suitable for an ILU coarse solve when using
-- AlgebraType("CPU", 1), so use the pivoting sparse direct SuperLU solver on
-- the coarse grid while retaining algebra type 1 on every level.
AssertPluginsLoaded({"NavierStokes", "SuperLU6"})

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3

-- Always refine the original UGX written by the mesh-generation script.
-- Opening and re-saving that file in ProMesh can preserve the serialized
-- NeuriteProjector while dropping its required npSurfParams attachment.
gridS ="third_trace_mm" 
gridName = "runs/"..gridS..".ugx"

numRefs = util.GetParamNumber("-numRefs", 0)
numPreRefs = 0 -- 0
v_max = 8.4 --mm/s
-- Physical parameters
viscosity 	= util.GetParamNumber("-visc", 1, "kinematic viscosity") -- kinematic viscosity [mm^2/s]
inflow		= util.GetParamNumber("-inflow", v_max, "max. inflow velocity")
inletRadius = util.GetParamNumber("-inletRadius", 1.0, "inlet lumen radius [mm]")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-noLaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactJac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type (no, full, weighted, lps, pos, reg)")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")


endTime = 0.02 -- [s]
dt = 0.01      -- [s]

local endTimeTag = tostring(endTime):gsub("%.", "p")
vtk_file_name = util.GetParam("-vtk",
    "Results/" .. gridS .. "_t" .. endTimeTag .. "_projector_ref" .. numRefs,
    "VTK output prefix")

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- Planned Darcy/BJS parameters. These are not assembled by the baseline.
-- The third_trace_mm mesh coordinates are in millimetres.
local permeabilityM2 = util.GetParamNumber("-permeabilityM2", 1.0e-12)
local permeabilityMm2 = permeabilityM2 * 1.0e6
local bjsAlpha = util.GetParamNumber("-bjsAlpha", 0.2)
local outerPressurePa = util.GetParamNumber("-outerPressurePa", 0.0)
assert(permeabilityM2 > 0 and bjsAlpha > 0)
-- D_na below is a legacy inactive value in um^2/s; convert when enabled.
D_na = 1.33e3

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

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

assert(file_contains(gridName, '<projector type="neurite"'),
    "The selected UGX has no serialized NeuriteProjector: " .. gridName)
assert(file_contains(gridName, 'name="npSurfParams"'),
    "The selected UGX has no npSurfParams attachment and cannot use projector " ..
    "refinement. Regenerate it with the mesh-generation Lua script and do not " ..
    "re-save it from ProMesh: " .. gridName)

print("Projector grid: " .. gridName)

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

-- The UGX contains detached surface faces. They do not participate in the
-- volume equation, so the unconnected-side check is skipped.
--[[
dom = util.CreateDomain(gridName,numPreRefs,
    {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
     "OuterWall", "Inlet", "Outlet"},
    true)
]]

-- This executable's neuro_collection plugin library is not present at runtime,
-- even though it is enabled in CMakeCache.txt.  On the first UGX load the core
-- NeuriteProjector constructor registers the npSurfParams attachment type, but
-- that happens after the grid reader has already skipped its stored values.
-- Preloading once registers the type so the real load below restores both the
-- attachment values and the serialized projector correctly.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {"Lumen", "Apical", "Inlet", "Outlet"}),
       "The domain is missing a required lumen subset")

print("Refinement mode: NEURITE PROJECTOR")
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
--approxSpace:add_fct("c", "Lagrange", 1)
--approxSpace:add_fct ({"u", "v","w", "p"}, "Lagrange", 1,"Lumen,Inlet,Outlet,Apical")
for _, fct in ipairs({"u", "v", "w", "p"}) do
    approxSpace:add_fct(fct, "Lagrange", 1, "Lumen,Inlet,Outlet,Apical")
end
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

-- Order the DoFs:
-- OrderLex crashes in DoFDistribution::permute_indices for this approximation
-- space in the selected UG4 build.  The default ordering is valid.
-- OrderLex(approxSpace, "xy")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------
--[[
dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c", "Apical")

diffusionEq = ConvectionDiffusionFV1("c", "Lumen,Membrane,Inter")
diffusionEq:set_diffusion(D_na)


]]--
--------------------------------------------------------------------------------
-- flow Equations 
--------------------------------------------------------------------------------

NavierStokesDisc = NavierStokesFV1 ({"u", "v","w", "p"}, {"Lumen"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (true)
NavierStokesDisc:set_laplace ( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity (viscosity)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

-- Unit tangent from the first two points of the centerline used to generate
-- this mesh; it points from the inlet into the lumen.
inlet_direction_x = -0.771159
inlet_direction_y = -0.349961
inlet_direction_z =  0.524920

function inletVelocity(x, y, z, t)
    return inflow * inlet_direction_x,
           inflow * inlet_direction_y,
           inflow * inlet_direction_z
end

InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inletVelocity", "Inlet")

OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Outlet")

WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Apical")

domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)

--------------------------------------------------------------------------------
-- MPI-friendly iterative solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace

-- The fine-grid iteration remains distributed.  SuperLU is used only for the
-- level-0 coarse problem, where scalar ILU fails on the velocity-pressure
-- saddle-point matrix.
solverDesc = {
    type = "newton",

    convCheck = {
        type = "standard",
        iterations = 100,
        absolute = 1e-15,
        reduction = 1e-10,
        verbose = true
    },

    linSolver = {
        type = "bicgstab",

        precond = {
            type = "gmg",
            smoother = {
                type = "gs", -- "jac",
                overlap = true,
                --damping = 0.66
            },
            cycle = "V",
            preSmooth = 4, -- 3,
            postSmooth = 4, -- 3,
            rap = true, -- false,
            baseLevel = numPreRefs,

            baseSolver = "superlu",
            gatheredBaseSolverIfAmbiguous = true
        },

        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-16,
            reduction = 0.5e-10, -- 1e-6,
            verbose = true
        }
    }
}

solver = util.solver.CreateSolver(solverDesc)

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection()
out:select_nodal ({"u", "v","w"}, "vel")
out:select_nodal ("u", "vel_u")
out:select_nodal ("v", "vel_v")
out:select_nodal ("w", "vel_w")
out:select_nodal ("p", "p")

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
--timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
for _, fct in ipairs({"u", "v", "w", "p"}) do
    Interpolate(0.0, u, fct)
end

--out:print(vtk_file_name .. "_init", u)

out:print_subsets(vtk_file_name .. "_lumen",u,"Lumen",0,0.0)
vtkObserver = LuaCallbackObserver()

function vtkCallback(step, time, currentDt)
    local currentSolution = vtkObserver:get_current_solution()

    out:print_subsets(
        vtk_file_name .. "_lumen",
        currentSolution,
        "Lumen",
        step,
        time
    )

    return 1
end

vtkObserver:set_callback("vtkCallback")
timeIntegrator:attach_observer(vtkObserver)


solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)

out:write_time_pvd(vtk_file_name .. "_lumen", u)

solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time/60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time/60.0))
