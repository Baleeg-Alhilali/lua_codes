-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

----------------------------------------------------------
--
--  One-concentration diffusion through the complete 3D geometry.
--
--  Grid: 1 um x 1 um rectangle.
--  Boundary conditions:
--    c1 = 1 on LeftBnd
--    c2 = 1 on RightBnd
--    top and bottom are natural no-flux walls
--    phi = 0 on LeftBnd as reference potential
--
----------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3
gridName = "runs/3rd_trace_243_third.ugx"

numRefs = 0
numPreRefs = 0


endTime = 0.5 -- [s]
dt = 0.005 
--dt = 1.0e-14      -- [s]

-- upwind = "partial"
vtk_file_name = "results/Comsol_MOdel" .. endTime.. numRefs

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

Faraday = 96485.33212       -- [C/mol]
gas_constant = 8.314462618  -- [J/(mol K)] check Joules
Temperature = 298.15        -- [K]

Porosity = 1.0
Viscosity = 1.0
Density = 1.0

D_na = 1.33e-9 -- Na+ diffusion coefficient in water [mm^2/s]
D_cl = 2.03e-9 -- Cl- diffusion coefficient in water [mm^2/s]

z1 = 1.0
z2 = -1.0

-- DarcyVelocityLinker uses K / mu. With mu = 1, K stores z F D / (R T).
MobilityNa = z1 * Faraday * D_na / (gas_constant * Temperature)
MobilityCl = z2 * Faraday * D_cl / (gas_constant * Temperature)

vacuumPermittivity = 8.8541878128e-15 -- [F/mm] --8.8541878128e-15 -- [F/mm]
relativePermittivity = 78.5
epsilon = relativePermittivity * vacuumPermittivity

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

-- The UGX contains six detached surface faces. They do not participate in the
-- volume diffusion equation, so skip the loader's unconnected-side check.
dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
     "OuterWall", "Inlet", "Outlet"},
    true)


balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

--OrderLex(approxSpace, "y")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------


dirichletBND = DirichletBoundary()
--dirichletBND:add(1.0e-9, "c", "Inlet")
--dirichletBND:add(0.0, "c", "Outlet")
dirichletBND:add(1.0, "c", "OuterWall")
--dirichletBND:add(1e-9, "c", "Baselateral")


EqS = 1.0e6
DiffusionEq = ConvectionDiffusionFV1("c", "Lumen,Membrane,Inter")
DiffusionEq:set_diffusion(D_na*EqS)
DiffusionEq:set_mass_scale(Porosity * EqS)


domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(DiffusionEq)
domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace

solverDesc = {
    type = "newton",

    lineSearch = {
        type = "standard",
        maxSteps = 8,
        lambdaStart = 1.0,
        lambdaReduce = 0.5,
        acceptBest = true,
        checkAll = false
    },

    convCheck = {
        type = "standard",
        iterations = 20,
        absolute = 1.0e-14,
        reduction = 1.0e-12,
        verbose = true
    },

	linSolver = {
    type = "lu",
    showProgress = true,
    info = true
}


}

solver = util.solver.CreateSolver(solverDesc)

--[[
util.solver.defaults.approxSpace = approxSpace
solverDesc = {
	type = "newton",

	lineSearch = {
		type = "standard",
		maxSteps = 8,
		lambdaStart = 1.0,
		lambdaReduce = 0.5,
		acceptBest = true,
		checkAll = false
	},

	convCheck = {
		type = "standard",
		iterations = 30,
		absolute = 1e-12,
		reduction = 1e-8,
		verbose = true
	},

	linSolver = {
		type = "bicgstab",
		precond = {
			type = "gmg",
			smoother = {
    			type = "ilu",
    			inversionEps = 1.0e-14,
   				 overlap = true,
   				 consistentInterfaces = true
				},
			cycle = "W",
			preSmooth = 2,
			postSmooth = 2,
			rap = false,
			baseLevel = numPreRefs,
			baseSolver = "lu"
		},

		convCheck = {
			type = "standard",
			iterations = 100,
			
			absolute = 1e-12,
			reduction = 1e-8, -- 1e-6,
			verbose = true
		}
	}
}

solver = util.solver.CreateSolver(solverDesc)
]]--
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


Interpolate(0.0,  u, "c")


out:print(vtk_file_name.."_init", u)

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)
solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f seconds", solve_elapsed_time))
print(string.format("Total script runtime:      %.6f seconds", total_elapsed_time))
