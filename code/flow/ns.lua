--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Solves the steady two-dimensional flow benchmark around a cylinder with a
-- stabilized vertex-centered Navier--Stokes discretization.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolves model-relative paths and loads the UG4 flow and balancing utilities.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")
AssertPluginsLoaded({"NavierStokes"})

------------------------------------------------------------------------------------------
-- Parameters
------------------------------------------------------------------------------------------
-- Only -numRefs and -output are command-line inputs. Edit all physical,
-- discretization, mesh, and solver parameters in this block.

-- Physical parameters
viscosity 	= 1e-3
inflow		= 0.5
bStokes 	= false

-- Numerical parameters of the discretization
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= 0
bNoLaplace 	= false
bExactJac 	= false
bPecletBlend= false
upwind      = "lps"
bPac        = false
stab        = "fields"
diffLength  = "cor"
outputFileName = util.GetParam(
	"-output", "Results/Flow-lev" .. numRefs,
	"VTK output prefix")

-- Parameters of the solver
dim 		= 2
gridName = "cylinderp.ugx"
inflow_y_0	= 0 -- y coord. of the lower corner of the inflow
inflow_y_1	= 0.41 -- y coord. of the upper corner of the inflow

-- Grid file name

print (" Geometry: " .. gridName .. ", dim = " .. dim)
print (" Physical parameter:")
print ("    visc		= " .. viscosity)
print ("    inflow		= " .. inflow)
print ("    Stokes		= " .. tostring (bStokes))
print (" Numerical parameter of the discretization:")
print ("    numRefs		= " .. numRefs)
print ("    numPreRefs	= " .. numPreRefs)
print ("	noLaplace	= " .. tostring (bNoLaplace))
print ("	exactJac	= " .. tostring (bExactJac))
print ("    PecletBlend = " .. tostring (bPecletBlend))
print ("    upwind		= " .. upwind)
print ("    pac			= " .. tostring (bPac))
print ("    stab		= " .. stab)
print ("    difflength	= " .. diffLength)
print (" Numerical parameter of the solvers:")

------------------------------------------------------------------------------------------
-- Mesh and approximation space
------------------------------------------------------------------------------------------
-- Loads/refines the benchmark mesh and defines first-order velocity/pressure.

InitUG (dim, AlgebraType("CPU", dim + 1))

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs, {"Inlet", "Outlet", "Inner", "UpperWall", "LowerWall", "CylinderWall"})
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)
approxSpace:add_fct ({"u", "v", "p"}, "Lagrange", 1)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace

------------------------------------------------------------------------------------------
-- Discretization and boundary conditions
------------------------------------------------------------------------------------------
-- Defines the parabolic inlet functions first and then assembles the interior,
-- inlet, outlet, and no-slip wall terms.

-- inner space
NavierStokesDisc = NavierStokesFV1 ({"u", "v", "p"}, {"Inner"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_laplace ( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity (viscosity)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

-- boundary condition at the inlet
inflow_R = (inflow_y_1 - inflow_y_0) / 2
function inflowVel_u (x, y, t)
	return inflow * (inflow_y_1 - y) * (y - inflow_y_0) / (inflow_R * inflow_R)
end
function inflowVel2d (x, y, t)
	return inflowVel_u (x, y, t), 0.0
end
InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inflowVel2d", "Inlet")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Outlet")

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("UpperWall,LowerWall,CylinderWall")

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)

-- create the assembled operator for the solver
assembledOp = AssembledOperator (domainDisc)

------------------------------------------------------------------------------------------
-- Solver
------------------------------------------------------------------------------------------
-- Uses Newton with a BiCGStab/GMG linear solve.

solverDesc =
{
	type = "newton", -- alternative for linear Stokes: a linear solver
	linSolver =
	{
		type = "bicgstab", -- alternative: "gmres"
		precond =
		{
			type = "gmg", -- alternative: "ilu" for a single-level solve
			rap = false,
			smoother =
			{
				type = "ilu",
					beta = -0.2, -- alternative: 0.0 for standard ILU
				consistentInterfaces = true
			},
			preSmooth = 1,
			postSmooth = 1,
				baseSolver = "lu", -- alternative: "superlu" when available
				baseLevel = 0
		},
		convCheck =
		{
			type		= "standard",
			iterations	= 128,
			absolute	= 1e-14,
			reduction	= 0.5e-8,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		= 5,
		lambdaStart		= 1,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= 256,
		absolute	= 1e-14,
		reduction	= 1e-8,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Time discretization and initial values
------------------------------------------------------------------------------------------
-- This is a steady model. The solution starts from zero with the analytical
-- inlet profile interpolated into the initial x-velocity.

-- grid function for the solution
u = GridFunction (approxSpace)

-- set the initial guess
u:set (0)
Interpolate (inflowVel_u, u, "u")

-- order the DoFs:
--OrderCuthillMcKee (approxSpace, true)
OrderLex (approxSpace,  "x")

-- initialize the solver
solver:init (assembledOp)
solver:prepare (u)

-- Apply the stationary solver.
if not solver:apply (u) then
	print ("The solver failed.")
	exit ()
end

------------------------------------------------------------------------------------------
-- Output
------------------------------------------------------------------------------------------
-- Writes velocity components and pressure to the requested VTK prefix.

out = VTKOutput()
out:clear_selection ()
out:select_nodal ({"u", "v"}, "vel")
out:select_nodal ("u", "vel_u")
out:select_nodal ("v", "vel_v")
out:select_nodal ("p", "p")

vtk_file_name = outputFileName
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end
print ("Output to file " .. vtk_file_name .. ".vtu")
out:print (vtk_file_name, u)

-- End of File
