-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Single-process fixed-flow urea transport using backward (implicit) Euler.
-- The stationary Stokes field is solved once before the transport time loop.

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

AssertPluginsLoaded({"NavierStokes", "SuperLU6", "neuro_collection"})

totalStartTime = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3
gridName = "runs/third_trace_mm.ugx"

numRefs = 1
numPreRefs = 0

endTime = 60.0 -- [s]
dt = 0.5       -- [s]

results_dir = "Results"
endTimeTag = tostring(endTime):gsub("%.", "p")
vtk_file_name = results_dir .. "/third_trace_mm_coupled_gmg_t" .. endTimeTag .. "_ref" .. numRefs

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- C_star = C / concentrationScale. The physical inlet remains 5e-6 mol/mm^3
-- while the numerical inlet value is exactly one.
concentrationScale = 5.0e-6
normalizedInletConcentration = 1.0

meanInletVelocity = 1.5            -- [mm/s]
peakInletVelocity = 2.0 * meanInletVelocity
inletRadius = 1.0e-3               -- [mm]
viscosity = 1.0                    -- [mm^2/s]

ureaDiffusion = 1.67e-3            -- [mm^2/s]
apicalPermeability = 3.78e-4       -- [mm/s]
basolateralPermeability = 3.70e-4  -- [mm/s]
membraneTransportDensity = 1.0

--------------------------------------------------------------------------------
-- Domain hierarchy and distribution
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
assert(NumProcs() == 1, "This implicit-Euler configuration must be run with exactly one MPI process")
print("MPI processes: " .. NumProcs())
print("Coarse solver grid: " .. gridName)
print("NeuriteProjector hierarchy refinements: " .. numRefs)
print(string.format("Concentration normalization: C_star = C / %.8g mol/mm^3", concentrationScale))

local gridFile = io.open(gridName, "r")
assert(gridFile ~= nil, "Cannot open grid file: " .. gridName)
gridFile:close()

-- Register the serialized neurite-projector attachment before the real load.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)

-- Older coarse UGX files mislabeled the two one-sided terminal annuli as
-- Basolateral. Repair them before refinement so Basolateral contains only the
-- two-sided Membrane--Inter interface and cannot create unsupported DoFs.
if ProcRank() == 0 then
    repair_nephron_transport_boundary_subsets(dom)
end

assert(util.CheckSubsets(dom, {"Lumen", "Membrane", "Basolateral", "Apical", "Inter", "OuterWall", "Inlet", "Outlet", "TripleJunction"}),
       "The grid is missing a required flow or transport subset")

-- Refine with the serialized NeuriteProjector. This keeps new Apical,
-- Basolateral, terminal-cap and TripleJunction vertices on the curved nephron
-- geometry while retaining the parent/child hierarchy required by GMG.
balancer.firstDistLvl = numRefs
balancer.Rebalance(dom)
if numRefs > numPreRefs then
    local projectorRefiner = GlobalDomainRefiner(dom)
    for level = numPreRefs + 1, numRefs do
        TerminateAbortedRun()
        projectorRefiner:refine()
        TerminateAbortedRun()

        -- Refinement creates edge-midpoint vertices.  Correct any inherited
        -- one-sided Basolateral labels before that level is distributed so a
        -- compartment never receives unsupported interface DoFs.
        if ProcRank() == 0 then
            repair_nephron_transport_boundary_subsets(dom)
        end
        balancer.Rebalance(dom)
    end
    delete(projectorRefiner)
end

print("Domain hierarchy:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- Steady Stokes block
--------------------------------------------------------------------------------

flowSpace = ApproximationSpace(dom)
for _, fct in ipairs({"u", "v", "w", "p"}) do
    flowSpace:add_fct(fct, "Lagrange", 1, "Lumen,Inlet,Outlet,Apical,TripleJunction")
end
flowSpace:init_levels()
flowSpace:init_top_surface()

flowDisc = NavierStokesFV1({"u", "v", "w", "p"}, {"Lumen"})
flowDisc:set_stokes(true)
flowDisc:set_laplace(true)
flowDisc:set_kinematic_viscosity(viscosity)
flowDisc:set_upwind("lps")
flowDisc:set_peclet_blend(false)
flowDisc:set_stabilization("fields", "cor")
flowDisc:set_pac_upwind(false)

local inletCenterX = 0.320613055
local inletCenterY = 0.471824404
local inletCenterZ = 0.242263399
local inletDirectionX = -0.7784261677084318
local inletDirectionY = -0.3468720130022301
local inletDirectionZ =  0.5231945221641233

function coupledGMGInletVelocity(x, y, z, t)
    local dx = x - inletCenterX
    local dy = y - inletCenterY
    local dz = z - inletCenterZ
    local axial = dx * inletDirectionX
                + dy * inletDirectionY
                + dz * inletDirectionZ
    local radialSquared = math.max(0.0, dx * dx + dy * dy + dz * dz - axial * axial)
    local relativeRadiusSquared = radialSquared / (inletRadius * inletRadius)
    local profile = math.max(0.0, 1.0 - relativeRadiusSquared)
    if relativeRadiusSquared >= 1.0 - 1.0e-6 then
        profile = 0.0
    end
    local speed = peakInletVelocity * profile
    return speed * inletDirectionX,
           speed * inletDirectionY,
           speed * inletDirectionZ
end

flowInlet = NavierStokesInflow(flowDisc)
flowInlet:add("coupledGMGInletVelocity", "Inlet")
flowOutlet = NavierStokesNoNormalStressOutflow(flowDisc)
flowOutlet:add("Outlet")
flowWall = NavierStokesWall(flowDisc)
flowWall:add("Apical,TripleJunction")

flowDomainDisc = DomainDiscretization(flowSpace)
flowDomainDisc:add(flowDisc)
flowDomainDisc:add(flowInlet)
flowDomainDisc:add(flowOutlet)
flowDomainDisc:add(flowWall)

local flowSmoother = ILU()
flowSmoother:set_damp(0.7)
flowSmoother:set_inversion_eps(1.0e-14)
flowSmoother:enable_overlap(true)
flowSmoother:enable_consistent_interfaces(true)

local flowGMG = GeometricMultiGrid(flowSpace)
flowGMG:set_base_level(numPreRefs)
flowGMG:set_base_solver(AgglomeratingSolver(SuperLU()))
flowGMG:set_gathered_base_solver_if_ambiguous(true)
flowGMG:set_smoother(flowSmoother)
flowGMG:set_cycle_type("V")
flowGMG:set_num_presmooth(4)
flowGMG:set_num_postsmooth(4)
flowGMG:set_rap(false)
flowGMG:set_transfer(StdTransfer())

util.solver.defaults.approxSpace = flowSpace
flowSolverDesc = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 50,
        absolute = 1.0e-14, reduction = 1.0e-12,
        verbose = true
    },
    linSolver = {
        type = "bicgstab", precond = flowGMG,
        convCheck = {
            type = "standard", iterations = 150,
            absolute = 1.0e-14, reduction = 1.0e-12,
            verbose = true
        }
    }
}

flowSolver = util.solver.CreateSolver(flowSolverDesc)

flowSolution = GridFunction(flowSpace)
for _, fct in ipairs({"u", "v", "w", "p"}) do
    Interpolate(0.0, flowSolution, fct)
end

-- Solve the stationary Stokes problem exactly once.  flowSolution is not
-- changed after this point; every transport step reads this same field.
print("Solving the stationary Stokes field once ...")
local flowSolveStartTime = GetClockS()
flowSolver:init(AssembledOperator(flowDomainDisc))
flowSolver:prepare(flowSolution)
assert(flowSolver:apply(flowSolution), "Stationary Stokes solve failed")
print(string.format("Stationary Stokes runtime: %.6f min", (GetClockS() - flowSolveStartTime) / 60.0))

-- The user data keeps a reference to flowSolution, so each new Stokes result
-- is immediately used by the transport discretization.
local solvedFlowVelocity = ExplicitGridFunctionVector(flowSolution, "u,v,w")

--------------------------------------------------------------------------------
-- Urea transport block
--------------------------------------------------------------------------------

transportSpace = ApproximationSpace(dom)
transportSpace:add_fct("urea_lumen", "Lagrange", 1, "Lumen,Apical,Inlet,Outlet,TripleJunction")
transportSpace:add_fct("urea_membrane", "Lagrange", 1, "Membrane,Apical,Basolateral,TripleJunction")
transportSpace:add_fct("urea_inter", "Lagrange", 1, "Inter,Basolateral,OuterWall,Inlet,Outlet,TripleJunction")
transportSpace:init_levels()
transportSpace:init_top_surface()

transportBoundary = DirichletBoundary()
transportBoundary:add(normalizedInletConcentration, "urea_lumen", "Inlet")
transportBoundary:add(0.0, "urea_lumen", "Outlet")
transportBoundary:add(0.0, "urea_inter", "OuterWall,Inlet,Outlet")

lumenTransport = ConvectionDiffusionFV1("urea_lumen", "Lumen")
lumenTransport:set_diffusion(ureaDiffusion)
lumenTransport:set_velocity(solvedFlowVelocity)
lumenTransport:set_upwind(FullUpwind())

membraneDiffusion = ConvectionDiffusionFV1("urea_membrane", "Membrane")
membraneDiffusion:set_diffusion(ureaDiffusion)
interDiffusion = ConvectionDiffusionFV1("urea_inter", "Inter")
interDiffusion:set_diffusion(ureaDiffusion)

apicalLeak = Leak({"urea_lumen", "urea_membrane"})
apicalLeak:set_permeability(apicalPermeability)
apicalTransport = MembraneTransportFV1("Apical", apicalLeak)
apicalTransport:set_density_function(membraneTransportDensity)

basolateralLeak = Leak({"urea_membrane", "urea_inter"})
basolateralLeak:set_permeability(basolateralPermeability)
basolateralTransport = MembraneTransportFV1("Basolateral", basolateralLeak)
basolateralTransport:set_density_function(membraneTransportDensity)

-- There is deliberately no urea_inter/Apical coupling.
transportDomainDisc = DomainDiscretization(transportSpace)
transportDomainDisc:add(lumenTransport)
transportDomainDisc:add(membraneDiffusion)
transportDomainDisc:add(interDiffusion)
transportDomainDisc:add(apicalTransport)
transportDomainDisc:add(basolateralTransport)
transportDomainDisc:add(transportBoundary)

util.solver.defaults.approxSpace = transportSpace
transportSolverDesc = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 100,
        absolute = 1.0e-15,
        reduction = 1.0e-15,
        verbose = true
    },
    linSolver = "superlu"
}

transportSolver = util.solver.CreateSolver(transportSolverDesc)

-- Theta = 1 gives first-order backward (implicit) Euler time discretization.
timeDisc = ThetaTimeStep(transportDomainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(transportSolver)
timeIntegrator:set_time_step(dt)

ureaSolution = GridFunction(transportSpace)
Interpolate(0.0, ureaSolution, "urea_lumen")
Interpolate(0.0, ureaSolution, "urea_membrane")
Interpolate(0.0, ureaSolution, "urea_inter")
-- Make the initial state consistent with the inlet/outlet Dirichlet data.
-- This prevents an artificial initial jump on constrained inlet DoFs.
transportDomainDisc:adjust_solution(ureaSolution)

--------------------------------------------------------------------------------
-- Urea-only ParaView output
--------------------------------------------------------------------------------

vtkLumen = VTKOutput()
vtkLumen:clear_selection()
vtkLumen:select("urea_lumen", "Urea_normalized")
vtkMembrane = VTKOutput()
vtkMembrane:clear_selection()
vtkMembrane:select("urea_membrane", "Urea_normalized")
vtkInter = VTKOutput()
vtkInter:clear_selection()
vtkInter:select("urea_inter", "Urea_normalized")

local outputTimes = {}
local lastOutputStep = 0

local function writeUreaVTK(step, time, solution)
    solution = solution or ureaSolution
    vtkLumen:print_subsets(vtk_file_name .. "_urea_lumen", solution, "Lumen", step, time)
    vtkMembrane:print_subsets(vtk_file_name .. "_urea_membrane", solution, "Membrane", step, time)
    vtkInter:print_subsets(vtk_file_name .. "_urea_inter", solution, "Inter", step, time)
    outputTimes[step] = time
    lastOutputStep = step
end

writeUreaVTK(0, 0.0)

--------------------------------------------------------------------------------
-- Fixed-flow transport with implicit Euler time stepping
--------------------------------------------------------------------------------

vtkObserver = LuaCallbackObserver()
function ureaVTKCallback(step, time, currentDt)
    writeUreaVTK(step, time, vtkObserver:get_current_solution())
    return 1
end
vtkObserver:set_callback("ureaVTKCallback")
timeIntegrator:attach_observer(vtkObserver)

transportSolver:init(AssembledOperator(transportDomainDisc))
transportSolver:prepare(ureaSolution)

local solveStartTime = GetClockS()
print(string.format("Implicit Euler transport: t=[0, %.8g] s with dt=%.8g s", endTime, dt))
assert(timeIntegrator:apply(ureaSolution, endTime, ureaSolution, 0.0), "Implicit Euler transport solve failed")

--------------------------------------------------------------------------------
-- PVD time series
--------------------------------------------------------------------------------

local function writeCompartmentPVD(prefix)
    if ProcRank() ~= 0 then return end
    local pvdFile = assert(io.open(prefix .. ".pvd", "w"))
    local basename = string.match(prefix, "([^/]+)$")
    local extension = NumProcs() > 1 and "pvtu" or "vtu"
    pvdFile:write('<?xml version="1.0"?>\n')
    pvdFile:write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n')
    for outputStep = 0, lastOutputStep do
        pvdFile:write(string.format('    <DataSet timestep="%.17g" group="" part="0" file="%s_t%04d.%s"/>\n',
                                    outputTimes[outputStep], basename, outputStep, extension))
    end
    pvdFile:write("  </Collection>\n</VTKFile>\n")
    pvdFile:close()
end

writeCompartmentPVD(vtk_file_name .. "_urea_lumen")
writeCompartmentPVD(vtk_file_name .. "_urea_membrane")
writeCompartmentPVD(vtk_file_name .. "_urea_inter")

local solveElapsedTime = GetClockS() - solveStartTime
local totalElapsedTime = GetClockS() - totalStartTime
print(string.format("Transport solve runtime: %.6f min", solveElapsedTime / 60.0))
print(string.format("Total runtime:         %.6f min", totalElapsedTime / 60.0))
print("ParaView array: Urea_normalized")
print(string.format("Physical urea [mol/mm^3] = Urea_normalized * %.8g", concentrationScale))
print("Urea ParaView output:")
print("  " .. vtk_file_name .. "_urea_lumen.pvd")
print("  " .. vtk_file_name .. "_urea_membrane.pvd")
print("  " .. vtk_file_name .. "_urea_inter.pvd")
