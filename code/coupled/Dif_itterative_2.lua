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

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3
gridS ="3rd_trace_third_boxref5" 
gridName = "ProMeshFiles/runs/"..gridS..".ugx"

numRefs = 2 --2
numPreRefs = 0 -- 0

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
dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
     "OuterWall", "Inlet", "Outlet"},
    true)

-- Refine and distribute the grid among all MPI processes.
-- use linear rifinment instead of the projector refinements
--local linearProjector = ProjectionHandler(dom:geometry3d(),dom:subset_handler())
-- dom:set_refinement_projector(linearProjector)

balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

-- Order the DoFs:
OrderLex (approxSpace, "xy")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c", "Apical")

diffusionEq = ConvectionDiffusionFV1("c", "Lumen,Membrane,Inter")
diffusionEq:set_diffusion(D_na)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(diffusionEq)
domainDisc:add(dirichletBND)

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
out:select_nodal("c", "concentration")

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
Interpolate(0.0, u, "c")

out:print(vtk_file_name .. "_init", u)

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)
solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time/60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time/60.0))
