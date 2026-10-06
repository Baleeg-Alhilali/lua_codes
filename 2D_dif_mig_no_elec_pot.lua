-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--------- Diffusion migration without an action potential unknown --- 
-- this model comes from the paper Benchmarks for multicomponent diffusion and electrochemical migration 
-- the idea is to link the migration term directly to gradient of the concentrations
-- --    d c_i / dt = - div(-D_i grad(c_i) + c_i v_i) + R_i


-------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------
dim = 2
gridName = "grids/charged_molecules_rect.ugx"
---
-- subset names ( Lumen , Memb , Inter , Blood,TopBottomWall,
--  Inter_blood,   LumenIn , BloodIn,Apical, Baseloteral)
numRefs = 3
numPreRefs = 1

endTime = 2.0e-3 -- [s]
dt = 1.0e-5      -- [s]

upwind = "partial"
results_dir = "results"
vtk_file_name = results_dir .. "/2D_simple_NaClK_" .. numRefs

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- Valences
zNa = 1.0
zCl = -1.0
zK = 1.0

-- Diffusion coefficients in water near room temperature [m^2/s].
D_na = 1.33e-9
D_cl = 2.03e-9
D_k = 1.96e-9

Faraday = 96485.33212       -- [C/mol]
gas_constant = 8.314462618  -- [J/(mol K)]
Temperature = 298.15        -- [K]

Porosity = 1.0
Viscosity = 1.0
Density = 1.0

Gravity = ConstUserVector(0.0)
--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 3))

dom = util.CreateDomain(gridName, numPreRefs, {"Inner", "LeftBnd", "BottomBnd", "TopBnd", "RightBnd"})
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())


approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("cNa", "Lagrange", 1)
approxSpace:add_fct("cCl", "Lagrange", 1)
approxSpace:add_fct("cK", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "y")

-----------------------------------
--- Setting up boundry conditions--
-----------------------------------
dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "cNa", "LeftBnd")
dirichletBND:add(1.0, "cCl", "RightBnd")
dirichletBND:add(0.5, "cK", "LeftBnd")


-------------------------------------
-- Setting the right hand side multiplier
-----------------------------------

function RHS(cNa,cCl,cK)
     return ((zNa * D_na * FlowEq1:gradient() + zCl * D_cl * FlowEq2:gradient() + zK * D_k * FlowEq3:gradient()))/(D_na*FlowEq1:value() + D_cl*FlowEq2:value() + D_k*FlowEq3:value())
end 
-----------------------------------
-- Setting up equations 
---------------------------------
FlowEq1 = ConvectionDiffusionFV1("cNa", "Inner")
FlowEq1:set_diffusion(D_na)
FlowEq1:set_source(zNa*FlowEq1:value()*RHS)

FlowEq2 = ConvectionDiffusionFV1("cCl", "Inner")
FlowEq2:set_diffusion(D_cl)
FlowEq2:set_source(zCl*cCl*RHS)

FlowEq3 = ConvectionDiffusionFV1("cK", "Inner")
FlowEq3:set_diffusion(D_k)
FlowEq3:set_source(zK*cK*RHS)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(FlowEq1)
domainDisc:add(FlowEq2)
domainDisc:add(FlowEq3)
domainDisc:add(dirichletBND)

-------------------------------------------------------------------------------
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
		iterations = 30,
		absolute = 1e-20,
		reduction = 1e-10,
		verbose = true
	},

	linSolver = {
		type = "bicgstab",
		precond = {
			type = "gmg",
			smoother = "ilu",
			cycle = "V",
			preSmooth = 2,
			postSmooth = 2,
			rap = false,
			baseLevel = numPreRefs,
			baseSolver = "lu"
		},

		convCheck = {
			type = "standard",
			iterations = 100,
			absolute = 1e-22,
			reduction = 1e-5,
			verbose = true
		}
	}
}

solver = util.solver.CreateSolver(solverDesc)

-------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection()
out:select_nodal("cNa", "cNa")
out:select_nodal("cCl", "cCl")
out:select_nodal("cK", "cK")


--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
Interpolate(0.0, u, "cNa")
Interpolate(0.0, u, "cCl")
Interpolate(0.0, u, "cK")

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

timeIntegrator:apply(u, endTime, u, 0.0)

