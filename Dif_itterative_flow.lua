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

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3
gridS ="third_trace_ogrid_no_projector" 
gridName = "runs/"..gridS..".ugx"

gridS ="4th_trace_third" 
gridName = "runs/"..gridS..".ugx"

v_max = 472e-7 --m/s
-- Physical parameters
viscosity 	= util.GetParamNumber("-visc", 1e-3, "kinematic viscosity")
inflow		= util.GetParamNumber("-inflow", v_max, "max. inflow velocity")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-noLaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactJac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type (no, full, weighted, lps, pos, reg)")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")


endTime = 2 -- [s]
dt = 0.1      -- [s]

vtk_file_name = "Results/"..gridS .. endTime .. "_ref" .. numRefs

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- The mesh coordinates are in micrometres, so diffusion is in um^2/s.
D_na = 1.33e3

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

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
NavierStokesDisc:set_stokes (bStokes)
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

-- This configuration deliberately contains no ILU, ILUT, or direct LU solve.
-- Damped Jacobi is inexpensive and parallel, while GMG controls the iteration
-- count as the mesh is refined.
solverDesc = {
    type = "newton",

    convCheck = {
        type = "standard",
        iterations = 5,
        absolute = 1e-8,
        reduction = 1e-8,
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

            -- Keep the coarse solve distributed and iterative. Using "lu"
            -- here would invoke the expensive ILUT factorization.
            baseSolver = {
                type = "bicgstab",
                precond = {
                    type = "ilu",
                    overlap = true,
                    --damping = 0.66
                },
                convCheck = {
                    type = "standard",
                    iterations = 200,
                    absolute = 1e-10,
                    reduction = 0.5e-8, -- 1e-4,
                    verbose = false
                }
            }
        },

        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 0.5e-8, -- 1e-6,
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
--[[
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
]]

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)

out:write_time_pvd(vtk_file_name .. "_lumen", u)

solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time/60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time/60.0))
