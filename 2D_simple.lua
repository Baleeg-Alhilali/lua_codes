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
--  2D simple coupled Na/Cl/K transport model.
--
--  Species equations are assembled with ConvectionDiffusionFV1:
--
--    d c_i / dt = - div(-D_i grad(c_i) + c_i v_i) + R_i
--
--  The migration part is moved to the right-hand side so the species equations
--  do not need DarcyVelocityLinker or a vector-valued Lua velocity callback.
--  In flux form:
--
--    Na: -D_Na(grad c_Na + c_Na G) + c_Na q
--    Cl: -D_Cl(grad c_Cl - c_Cl G) + c_Cl q
--    K : -D_K (grad c_K  + c_K  G) + c_K  q
--
--  is written as ordinary diffusion/advection plus a right-hand-side term:
--
--    d c_i / dt = -div(-D_i grad(c_i) + c_i q)
--                 + div(D_i z_i c_i G) + R_i
--
--  The algebraic electroneutral field used in that RHS term is
--
--    G = (q(c_Na - c_Cl + c_K)
--         - (D_Na grad c_Na - D_Cl grad c_Cl + D_K grad c_K))
--        / (D_Na c_Na + D_Cl c_Cl + D_K c_K)
--
----------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 2
gridName = "grids/charged_molecules_rect.ugx"

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
DNa = 1.33e-9
DCl = 2.03e-9
DK = 1.96e-9

-- Constant background advection q [m/s].  With q = 0 the remaining coupling is
-- the diffusive charge-compensation term in G.
qX = 0.0
qY = 0.0

-- Reactions. Keep them explicit so they can be replaced by uptake/secretion terms.
RNa = 0.0
RCl = 0.0
RK = 0.0

Porosity = 1.0
minDenominator = 1.0e-30

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 3))

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------
-- Helper functions for the coupled algebraic right-hand side
--------------------------------------------------------------------------------

function VecComp(v, i)
	if type(v) == "table" then
		return v[i] or v[i - 1] or 0.0
	end

	if v ~= nil then
		local okOneBased, oneBased = pcall(function()
			return v[i]
		end)
		if okOneBased and oneBased ~= nil then
			return oneBased
		end

		local okZeroBased, zeroBased = pcall(function()
			return v[i - 1]
		end)
		if okZeroBased and zeroBased ~= nil then
			return zeroBased
		end
	end

	return 0.0
end

function GComponents(cNa, cCl, cK, gradNa, gradCl, gradK)
	local denom = DNa * cNa + DCl * cCl + DK * cK
	if math.abs(denom) < minDenominator then
		denom = minDenominator
	end

	local charge = cNa - cCl + cK

	local dNaX = VecComp(gradNa, 1)
	local dNaY = VecComp(gradNa, 2)
	local dClX = VecComp(gradCl, 1)
	local dClY = VecComp(gradCl, 2)
	local dKX = VecComp(gradK, 1)
	local dKY = VecComp(gradK, 2)

	local gx = (qX * charge - (DNa * dNaX - DCl * dClX + DK * dKX)) / denom
	local gy = (qY * charge - (DNa * dNaY - DCl * dClY + DK * dKY)) / denom

	return gx, gy
end

function MigrationRhs(c, z, D, cNa, cCl, cK, gradNa, gradCl, gradK)
	local gx, gy = GComponents(cNa, cCl, cK, gradNa, gradCl, gradK)

	-- Conservative form is div(D z c G).  ConvectionDiffusionFV1 exposes only a
	-- scalar source here, so this keeps the moved term as a local RHS correction.
	-- For a strict finite-volume flux, UG4 needs a vector flux/linker object.
	return D * z * c * (gx + gy)
end

function SourceNaFct(cNa, cCl, cK, gradNa, gradCl, gradK)
	return RNa + MigrationRhs(cNa, zNa, DNa, cNa, cCl, cK, gradNa, gradCl, gradK)
end

function SourceClFct(cNa, cCl, cK, gradNa, gradCl, gradK)
	return RCl + MigrationRhs(cCl, zCl, DCl, cNa, cCl, cK, gradNa, gradCl, gradK)
end

function SourceKFct(cNa, cCl, cK, gradNa, gradCl, gradK)
	return RK + MigrationRhs(cK, zK, DK, cNa, cCl, cK, gradNa, gradCl, gradK)
end

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "cNa", "LeftBnd")
dirichletBND:add(1.0, "cCl", "RightBnd")
dirichletBND:add(1.0, "cK", "LeftBnd")

FlowNa = ConvectionDiffusionFV1("cNa", "Inner")
FlowNa:set_upwind(UpwindFV1(upwind))
FlowNa:set_mass_scale(Porosity)
FlowNa:set_diffusion(DNa)

FlowCl = ConvectionDiffusionFV1("cCl", "Inner")
FlowCl:set_upwind(UpwindFV1(upwind))
FlowCl:set_mass_scale(Porosity)
FlowCl:set_diffusion(DCl)

FlowK = ConvectionDiffusionFV1("cK", "Inner")
FlowK:set_upwind(UpwindFV1(upwind))
FlowK:set_mass_scale(Porosity)
FlowK:set_diffusion(DK)

SourceNa = LuaUserFunctionNumber("SourceNaFct", 6)
SourceCl = LuaUserFunctionNumber("SourceClFct", 6)
SourceK = LuaUserFunctionNumber("SourceKFct", 6)

SourceNa:set_input(0, FlowNa:value())
SourceNa:set_input(1, FlowCl:value())
SourceNa:set_input(2, FlowK:value())
SourceNa:set_input(3, FlowNa:gradient())
SourceNa:set_input(4, FlowCl:gradient())
SourceNa:set_input(5, FlowK:gradient())

SourceCl:set_input(0, FlowNa:value())
SourceCl:set_input(1, FlowCl:value())
SourceCl:set_input(2, FlowK:value())
SourceCl:set_input(3, FlowNa:gradient())
SourceCl:set_input(4, FlowCl:gradient())
SourceCl:set_input(5, FlowK:gradient())

SourceK:set_input(0, FlowNa:value())
SourceK:set_input(1, FlowCl:value())
SourceK:set_input(2, FlowK:value())
SourceK:set_input(3, FlowNa:gradient())
SourceK:set_input(4, FlowCl:gradient())
SourceK:set_input(5, FlowK:gradient())

FlowNa:set_source(SourceNa)
FlowCl:set_source(SourceCl)
FlowK:set_source(SourceK)

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(FlowNa)
domainDisc:add(FlowCl)
domainDisc:add(FlowK)
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

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

os.execute("mkdir -p " .. results_dir)

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
