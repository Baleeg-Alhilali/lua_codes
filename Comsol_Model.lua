-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
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

numRefs = 5
numPreRefs = 1


endTime = 6-- [s]
dt = 0.1
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

D_na = 1.33e-3 -- Na+ diffusion coefficient in water [mm^2/s]
D_cl = 2.03e-3 -- Cl- diffusion coefficient in water [mm^2/s]

z1 = 1.0
z2 = -1.0

-- DarcyVelocityLinker uses K / mu. With mu = 1, K stores z F D / (R T).
MobilityNa = z1 * Faraday * D_na / (gas_constant * Temperature)
MobilityCl = z2 * Faraday * D_cl / (gas_constant * Temperature)

vacuumPermittivity = 8.8541878128e-12 -- [F/mm] --8.8541878128e-15 -- [F/mm]
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
dirichletBND:add(0.0, "c2", "LeftBnd,TopBnd,BottomBnd")
dirichletBND:add(1.0, "c2", "RightBnd")
dirichletBND:add(0.0, "c1", "RightBnd,TopBnd,BottomBnd")
dirichletBND:add(0.0, "phi", "TopCorner")
--dirichletBND:add(z2 *Faraday, "phi", "RightBnd")

--[[

PotentialEq = ConvectionDiffusionFV1("phi", "Inner")
FlowEq1 = ConvectionDiffusionFV1("c1", "Inner")
FlowEq1:set_upwind(FullUpwind()) -- (UpwindFV1(upwind))
FlowEq1:set_mass_scale(Porosity)
FlowEq1:set_flux(z1*D_na*Faraday/(gas_constant*Temperature) * FlowEq1:value()*PotentialEq:gradient())
FlowEq1:set_diffusion(D_na)

FlowEq2 = ConvectionDiffusionFV1("c2", "Inner")
FlowEq2:set_upwind(FullUpwind()) -- (UpwindFV1(upwind))FullUpwind()
FlowEq2:set_mass_scale(Porosity)
FlowEq2:set_flux(z2*D_cl*Faraday/(gas_constant*Temperature) * FlowEq2:value()*PotentialEq:gradient())
FlowEq2:set_diffusion(D_cl)



PotentialEq:set_mass_scale(0.0)
PotentialEq:set_diffusion((Faraday * Faraday)/(gas_constant * Temperature) *( (z1 * z1 * D_na * FlowEq1:value()) + (z2 * z2 * D_cl * FlowEq2:value())))
PotentialEq:set_flux(Faraday * ( (z1 * D_na * FlowEq1:gradient()) + (z2  * D_cl * FlowEq2:gradient())))

]]

PotentialEq = ConvectionDiffusionFV1("phi", "Inner")
FlowEq1 = ConvectionDiffusionFV1("c1", "Inner")
FlowEq2 = ConvectionDiffusionFV1("c2", "Inner")

PotentialEq:set_mass_scale(0.0)
PotentialEq:set_diffusion(0.0)

NegativeConductivity = ScaleAddLinkerNumber()
NegativeConductivity:add(
    (Faraday * Faraday * z1 * z1 * D_na)
        / (gas_constant * Temperature),
    FlowEq1:value()
)
NegativeConductivity:add(
    (Faraday * Faraday * z2 * z2 * D_cl)
        / (gas_constant * Temperature),
    FlowEq2:value()
)


FlowEq1:set_upwind(FullUpwind())
FlowEq1:set_mass_scale(Porosity)
FlowEq1:set_flux(
    -z1*D_na*Faraday/(gas_constant*Temperature)
    * FlowEq1:value()*PotentialEq:gradient()
)
FlowEq1:set_diffusion(D_na)

FlowEq2:set_upwind(FullUpwind())
FlowEq2:set_mass_scale(Porosity)
FlowEq2:set_flux(
    -z2*D_cl*Faraday/(gas_constant*Temperature)
    * FlowEq2:value()*PotentialEq:gradient()
)
FlowEq2:set_diffusion(D_cl)



PotentialFlux = ScaleAddLinkerVector()
PotentialFlux:add(NegativeConductivity, PotentialEq:gradient())
PotentialFlux:add(Faraday * z1 * D_na, FlowEq1:gradient())
PotentialFlux:add(Faraday * z2 * D_cl, FlowEq2:gradient())

PotentialEq:set_flux(PotentialFlux)

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
		acceptBest = false,
		checkAll = false
	},

	convCheck = {
		type = "standard",
		iterations = 30,
		absolute = 1e-13,
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
			
			absolute = 1e-13,
			reduction = 1e-6, -- 1e-6,
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
--out:select(MigrationVelocity1, "migration_velocity_c1")
--out:select(MigrationVelocity2, "migration_velocity_c2")

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name, out))
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
-- The conductivity in the phi equation is proportional to
-- D_na*c1 + D_cl*c2.  A zero concentration initial guess makes the potential
-- operator singular during the stationary initialization below.  Start from
-- a positive, electroneutral concentration instead; the Dirichlet conditions
-- still overwrite the prescribed boundary values.
Interpolate(0.10, u, "c1")
Interpolate(0.10, u, "c2")
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
