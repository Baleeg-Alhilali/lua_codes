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
--  Fixed charged molecule model:
--  Na+ and Cl- diffusion with electric migration through phi.
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

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 2
gridName = "grids/charged_molecules_rect.ugx"

numRefs = 6
numPreRefs = 1


endTime = 2.0e-5 -- [s]
dt = 1.0e-7
--dt = 1.0e-14      -- [s]

-- upwind = "partial"
vtk_file_name = "results/Charged_Mole_Fixed" .. numRefs

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

Faraday = 96485.33212       -- [C/mol]
gas_constant = 8.314462618  -- [J/(mol K)] check Joules
Temperature = 298.15        -- [K]

Porosity = 1.0
Viscosity = 1.0
Density = 1.0

D_na = 1.33e-3 -- Na+ diffusion coefficient in water [mm^2/s]
D_cl = 2.03e-3 -- Cl- diffusion coefficient in water [mm^2/s]

z1 = -1.0
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

InitUG(dim, AlgebraType("CPU", 3))

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = util.CreateDomain(gridName, numPreRefs, {"Inner", "LeftBnd", "BottomBnd", "TopBnd", "RightBnd","TopCorner"})
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("c1", "Lagrange", 1)
approxSpace:add_fct("c2", "Lagrange", 1)
approxSpace:add_fct("phi", "Lagrange", 1)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "y")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

Gravity = ConstUserVector(0.0)

dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c1", "LeftBnd")
dirichletBND:add(0.0, "c2", "LeftBnd,TopBnd,BottomBnd,TopCorner")
dirichletBND:add(1.0, "c2", "RightBnd")
dirichletBND:add(0.0, "c1", "RightBnd,TopBnd,BottomBnd,TopCorner")
dirichletBND:add(0.0, "phi", "TopCorner")
--dirichletBND:add(z2 *Faraday, "phi", "RightBnd")


---- Added fix by the agent
function ChargeDensityFct(c1, c2)
	return Faraday *1e-9 *(z1 * c1 + z2 * c2) -- 1e-9 to convert from mol/m^3 to mol/mm^3
end

function DChargeDensityFct_c1(c1, c2)
	return Faraday * z1 * 1e-9
end

function DChargeDensityFct_c2(c1, c2)
	return Faraday * z2 * 1e-9
end

function sigma(c1, c2)
	return Faraday / (R*T) * ((z1*D_na*c1)+(z2*D_cl*c2))
end

function RHS(c1_grad, c2_grad)
	return Faraday *  ((z1*D_na*c1_grad)+(z2*D_cl*c2_grad))
end


ChargeDensity = LuaUserFunctionNumber("ChargeDensityFct", 2)
ChargeDensity:set_deriv(0, "DChargeDensityFct_c1")
ChargeDensity:set_deriv(1, "DChargeDensityFct_c2")

ChargeDensity = LuaUserFunctionNumber("ChargeDensityFct", 2)
ChargeDensity:set_deriv(0, "DChargeDensityFct_c1")
ChargeDensity:set_deriv(1, "DChargeDensityFct_c2")
--
PotentialEq = ConvectionDiffusionFV1("phi", "Inner")
PotentialEq:set_mass_scale(0.0)
PotentialEq:set_diffusion(epsilon)
PotentialEq:set_source(ChargeDensity)

MigrationVelocity1 = DarcyVelocityLinker()
MigrationVelocity1:set_permeability(MobilityNa)
MigrationVelocity1:set_viscosity(Viscosity)
MigrationVelocity1:set_density(Density)
MigrationVelocity1:set_gravity(Gravity)

MigrationVelocity2 = DarcyVelocityLinker()
MigrationVelocity2:set_permeability(MobilityCl)
MigrationVelocity2:set_viscosity(Viscosity)
MigrationVelocity2:set_density(Density)
MigrationVelocity2:set_gravity(Gravity)

FlowEq1 = ConvectionDiffusionFV1("c1", "Inner")
FlowEq1:set_upwind(FullUpwind()) -- (UpwindFV1(upwind))
FlowEq1:set_mass_scale(Porosity)
FlowEq1:set_velocity(MigrationVelocity1)
FlowEq1:set_diffusion(D_na)

FlowEq2 = ConvectionDiffusionFV1("c2", "Inner")
FlowEq2:set_upwind(FullUpwind()) -- (UpwindFV1(upwind))FullUpwind()
FlowEq2:set_mass_scale(Porosity)
FlowEq2:set_velocity(MigrationVelocity2)
FlowEq2:set_diffusion(D_cl)

ChargeDensity:set_input(0, FlowEq1:value())
ChargeDensity:set_input(1, FlowEq2:value())
MigrationVelocity1:set_pressure_gradient(PotentialEq:gradient())
MigrationVelocity2:set_pressure_gradient(PotentialEq:gradient())

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(PotentialEq)
domainDisc:add(FlowEq1)
domainDisc:add(FlowEq2)
domainDisc:add(dirichletBND)

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
		iterations = 30,
		absolute = 1e-8,
		reduction = 1e-6,
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
			
			absolute = 1e-8,
			reduction = 1e-8, -- 1e-6,
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
out:select_nodal("c1", "c1")
out:select_nodal("c2", "c2")
out:select_nodal("phi", "phi")
out:select(MigrationVelocity1, "migration_velocity_c1")
out:select(MigrationVelocity2, "migration_velocity_c2")

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
Interpolate(0.0, u, "c1")
Interpolate(0.0, u, "c2")
Interpolate(0.0, u, "phi")

fixer = DirichletBoundary()
fixer:invert_subset_selection()
fixer:add("c1", "")
fixer:add("c2", "")
domainDisc:add(fixer)

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)
solver:apply(u)

domainDisc:remove(fixer)

out:print(vtk_file_name.."_init", u)

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

timeIntegrator:apply(u, endTime, u, 0.0)
