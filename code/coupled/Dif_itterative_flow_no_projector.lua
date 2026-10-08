-- Resolve every relative input and output path from the Model directory.
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
gridS ="3rd_trace_third_dani" --"third_trace_third" --"third_trace_third_dani"
gridName = "runs/"..gridS..".ugx"

numRefs = util.GetParamNumber("-numRefs", 0)
numPreRefs = 0 -- 0
v_max = 8.4e3 --um/s
-- Physical parameters
viscosity 	= util.GetParamNumber("-visc", 1e6, "kinematic viscosity")
inflow		= util.GetParamNumber("-inflow", v_max, "max. inflow velocity")
inletRadius = util.GetParamNumber("-inletRadius", 1.0, "inlet lumen radius [um]")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")
bNoLaplace 	= util.HasParamOption("-noLaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactJac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "lps", "Upwind type (no, full, weighted, lps, pos, reg)")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")


endTime = 0.2 -- [s]
dt = 0.1      -- [s]

vtk_file_name = "Results/"..gridS  .. "_ref" .. numRefs

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

-- This comparison mesh intentionally contains no serialized NeuriteProjector
-- and no npSurfParams attachment. Load it once, then explicitly install a plain
-- ProjectionHandler. Its default RefinementProjector performs midpoint/linear
-- refinement only and cannot move new vertices back onto the SWC tube.
dom = Domain()
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {"Lumen", "Apical", "Inlet", "Outlet"}),
       "The domain is missing a required lumen subset")

local linearProjector = ProjectionHandler(
    dom:geometry3d(), dom:subset_handler())
dom:set_refinement_projector(linearProjector)

print("Refinement mode: LINEAR (NeuriteProjector disabled)")
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
--OrderLex(approxSpace, "xy")
--OrderCuthillMcKee(approxSpace, true)

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

-- Fully developed circular Poiseuille profile at the inlet:
--
--     u(r) = vmax * (1 - r^2/R^2) * inlet_direction
--
-- Thus the centerline speed is vmax (= -inflow), the wall speed is zero, and
-- the analytical area-average speed is vmax/2.  The center and direction are
-- taken from the first inlet cross-section of third_trace_ogrid.  The axial
-- component is removed before evaluating r, so the profile remains robust to
-- tiny round-off deviations of cap vertices from the inlet plane.
inlet_center_x = 512.472
inlet_center_y = 427.021
inlet_center_z = 164.835

inlet_direction_x = -0.785486660321187
inlet_direction_y = -0.339816285667396
inlet_direction_z =  0.517238434817718

function inletVelocity(x, y, z, t)
    local dx = x - inlet_center_x
    local dy = y - inlet_center_y
    local dz = z - inlet_center_z

    local axial = dx * inlet_direction_x
                + dy * inlet_direction_y
                + dz * inlet_direction_z
    local radialSquared = dx * dx + dy * dy + dz * dz - axial * axial
    radialSquared = math.max(0.0, radialSquared)

    local relativeRadiusSquared = radialSquared / (inletRadius * inletRadius)
    local profile = math.max(0.0, 1.0 - relativeRadiusSquared)

    -- Make the shared Inlet/Apical rim exactly no-slip despite coordinate
    -- round-off in the serialized mesh.
    if relativeRadiusSquared >= 1.0 - 1.0e-6 then
        profile = 0.0
    end

    local speed = inflow * profile
    return speed * inlet_direction_x,
           speed * inlet_direction_y,
           speed * inlet_direction_z
end

print(string.format(
    "Inlet profile: Poiseuille, vmax=%.6g um/s, average=%.6g um/s, R=%.6g um",
    inflow, 0.5 * inflow, inletRadius))

InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inletVelocity", "Inlet")

OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Outlet")

WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Apical")

--WallDisc = DirichletBoundary()
--WallDisc:add (0, "u", "Apical")
--WallDisc:add (0, "v", "Apical")
--WallDisc:add (0, "w", "Apical")

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
        iterations = 5,
        absolute = 1e-2,
        reduction = 1e-8,
        verbose = true
    },

    linSolver = {
        type = "bicgstab",

        precond = {
            type = "gmg",
            smoother = {
                type = "ilu", -- "jac",
                overlap = true,
                --consistentInterfaces = true,
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
            reduction = 0.5e-16, -- 1e-6,
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

out:print(vtk_file_name .. "_init", u)

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
