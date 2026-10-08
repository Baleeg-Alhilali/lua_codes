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

AssertPluginsLoaded({"NavierStokes", "SuperLU6"})

local totalStartTime = GetClockS()

--------------------------------------------------------------------------------
-- Run parameters
--------------------------------------------------------------------------------

dim = 3

-- The concentration unknowns are normalized with C_ref so that their values
-- remain O(1).  The physical concentration is C = C_star * C_ref.
local concentrationScale = 1e-6 -- mol/mm^3
local normalizedInletConcentration = 5.0

local refinementLevel = util.GetParamNumber("-numRefs", 1)
assert(refinementLevel >= 0 and refinementLevel == math.floor(refinementLevel),
    "-numRefs must be a non-negative integer")

-- Refinement is deliberately not performed inside this script.  Direct UG4
-- refinement creates coincident but disconnected vertices at compartment
-- interfaces.  Those vertices produce unsupported concentration DoFs and NaN
-- defects.  A conforming, repaired mesh is therefore loaded for level > 0.
local defaultGridName
if refinementLevel == 0 then
    defaultGridName = "runs/third_trace_mm.ugx"
else
    defaultGridName = string.format(
        "runs/third_trace_mm_refined_%d_conforming.ugx",
        refinementLevel
    )
end

gridName = util.GetParam("-grid", defaultGridName, "input UGX grid")

local gridFile = io.open(gridName, "r")
assert(gridFile ~= nil,
    "Cannot open grid file: " .. gridName ..
    ". For refined runs, generate the conforming repaired mesh first.")
gridFile:close()

local endTime = util.GetParamNumber("-endTime", 60.0, "end time [s]")
local dt = util.GetParamNumber("-dt", 0.5, "time-step size [s]")
assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0, "-dt must be positive")

local meanInletVelocity = util.GetParamNumber(
    "-vavg", 1.5, "mean inlet velocity [mm/s]"
)
local peakInletVelocity = util.GetParamNumber(
    "-inflow", 2.0 * meanInletVelocity,
    "peak parabolic inlet velocity [mm/s]"
)
local inletRadius = util.GetParamNumber(
    "-inletRadius", 0.001, "lumen inlet radius [mm]"
)
local viscosity = util.GetParamNumber(
    "-visc", 1.0, "kinematic viscosity [mm^2/s]"
)
local flowUpwind = util.GetParam(
    "-upwind", "lps", "flow upwind type"
)
local flowStabilization = util.GetParam(
    "-stab", "fields", "flow stabilization type"
)
local flowDiffusionLength = util.GetParam(
    "-difflength", "cor", "flow diffusion-length type"
)

local endTimeTag = tostring(endTime):gsub("%.", "p")
local vtkPrefix = util.GetParam(
    "-vtk",
    string.format(
        "Results/third_trace_mm_t%s_ref%d_corrected",
        endTimeTag,
        refinementLevel
    ),
    "VTK output prefix"
)

--------------------------------------------------------------------------------
-- Physical transport parameters (mesh unit: mm)
--------------------------------------------------------------------------------

local ureaDiffusion = 1.67e-3             -- mm^2/s
local apicalPermeability = 3.78e-4        -- mm/s
local basolateralPermeability = 3.70e-4   -- mm/s
local membraneTransportDensity = 1.0

--------------------------------------------------------------------------------
-- Domain
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print("Grid: " .. gridName)
print("Requested pre-refined level: " .. refinementLevel)
print(string.format(
    "Concentration normalization: C_star = C / %.8g mol/mm^3",
    concentrationScale
))

-- Preloading registers the neurite-projector attachment types before the real
-- load.  It is harmless for a flat repaired grid and required by some builds.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)

assert(util.CheckSubsets(dom, {
    "Lumen", "Membrane", "Basolateral", "Apical", "Inter",
    "OuterWall", "Inlet", "Outlet", "TripleJunction"
}), "The domain is missing a required flow or transport subset")

-- Do not call RefineAndRebalanceDomain here.  The selected refined mesh is
-- already refined and repaired so that compartment interfaces are conforming.
print("Domain info:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- 1. Solve the steady Stokes velocity field
--------------------------------------------------------------------------------

flowSpace = ApproximationSpace(dom)
for _, fct in ipairs({"u", "v", "w", "p"}) do
    flowSpace:add_fct(
        fct, "Lagrange", 1,
        "Lumen,Inlet,Outlet,Apical,TripleJunction"
    )
end
flowSpace:init_levels()
flowSpace:init_top_surface()

flowDisc = NavierStokesFV1({"u", "v", "w", "p"}, {"Lumen"})
flowDisc:set_stokes(true)
flowDisc:set_laplace(true)
flowDisc:set_kinematic_viscosity(viscosity)
flowDisc:set_upwind(flowUpwind)
flowDisc:set_peclet_blend(false)
flowDisc:set_stabilization(flowStabilization, flowDiffusionLength)
flowDisc:set_pac_upwind(false)

-- Center of the inlet cross-section and inward unit tangent from the SWC.
local inletCenterX = 0.32061305499999998
local inletCenterY = 0.47182440400000003
local inletCenterZ = 0.24226339899999999
local inletDirectionX = -0.771159
local inletDirectionY = -0.349961
local inletDirectionZ =  0.524920

function correctedInletVelocity(x, y, z, t)
    local dx = x - inletCenterX
    local dy = y - inletCenterY
    local dz = z - inletCenterZ
    local axialDistance = dx * inletDirectionX
                        + dy * inletDirectionY
                        + dz * inletDirectionZ
    local radialSquared = math.max(
        0.0,
        dx * dx + dy * dy + dz * dz - axialDistance * axialDistance
    )
    local profile = math.max(
        0.0,
        1.0 - radialSquared / (inletRadius * inletRadius)
    )
    local speed = peakInletVelocity * profile
    return speed * inletDirectionX,
           speed * inletDirectionY,
           speed * inletDirectionZ
end

flowInlet = NavierStokesInflow(flowDisc)
flowInlet:add("correctedInletVelocity", "Inlet")

flowOutlet = NavierStokesNoNormalStressOutflow(flowDisc)
flowOutlet:add("Outlet")

flowWall = NavierStokesWall(flowDisc)
-- The triple-junction surface belongs to the no-slip lumen wall.
flowWall:add("Apical,TripleJunction")

flowDomainDisc = DomainDiscretization(flowSpace)
flowDomainDisc:add(flowDisc)
flowDomainDisc:add(flowInlet)
flowDomainDisc:add(flowOutlet)
flowDomainDisc:add(flowWall)

util.solver.defaults.approxSpace = flowSpace
print("Solving steady Stokes velocity field ...")
local _, flowSolution = util.solver.SolveLinearProblem(
    flowDomainDisc,
    "superlu"
)

-- Bridge the separately solved velocity into the transport discretization.
-- ConvectionDiffusionFV1 needs one vector-valued user-data object; passing the
-- Lua table {"u", "v", "w"} is not a valid set_velocity overload.
local solvedFlowVelocity = ExplicitGridFunctionVector(
    flowSolution,
    "u,v,w"
)

--------------------------------------------------------------------------------
-- 2. Solve normalized urea transport
--------------------------------------------------------------------------------

transportSpace = ApproximationSpace(dom)
transportSpace:add_fct(
    "urea_lumen", "Lagrange", 1,
    "Lumen,Apical,Inlet,Outlet,TripleJunction"
)
transportSpace:add_fct(
    "urea_membrane", "Lagrange", 1,
    "Membrane,Apical,Basolateral,TripleJunction"
)
transportSpace:add_fct(
    "urea_inter", "Lagrange", 1,
    "Inter,Basolateral,OuterWall,Inlet,Outlet,TripleJunction"
)
transportSpace:init_levels()
transportSpace:init_top_surface()

print("Transport approximation space:")
transportSpace:print_statistic()

transportDomainDisc = DomainDiscretization(transportSpace)

transportBoundary = DirichletBoundary()
transportBoundary:add(
    normalizedInletConcentration,
    "urea_lumen",
    "Inlet"
)
transportBoundary:add(0.0, "urea_lumen", "Outlet")
transportBoundary:add(0.0, "urea_inter", "OuterWall")

lumenTransport = ConvectionDiffusionFV1("urea_lumen", "Lumen")
lumenTransport:set_diffusion(ureaDiffusion)
lumenTransport:set_velocity(solvedFlowVelocity)
-- Full upwind prevents the convection-dominated solution from developing the
-- negative oscillations seen with an unstabilized central flux.
lumenTransport:set_upwind(FullUpwind())
transportDomainDisc:add(lumenTransport)

membraneDiffusion = ConvectionDiffusionFV1("urea_membrane", "Membrane")
membraneDiffusion:set_diffusion(ureaDiffusion)
transportDomainDisc:add(membraneDiffusion)

interDiffusion = ConvectionDiffusionFV1("urea_inter", "Inter")
interDiffusion:set_diffusion(ureaDiffusion)
transportDomainDisc:add(interDiffusion)

apicalLeak = Leak({"urea_lumen", "urea_membrane"})
apicalLeak:set_permeability(apicalPermeability)
apicalTransport = MembraneTransportFV1("Apical", apicalLeak)
apicalTransport:set_density_function(membraneTransportDensity)
transportDomainDisc:add(apicalTransport)

basolateralLeak = Leak({"urea_membrane", "urea_inter"})
basolateralLeak:set_permeability(basolateralPermeability)
basolateralTransport = MembraneTransportFV1(
    "Basolateral",
    basolateralLeak
)
basolateralTransport:set_density_function(membraneTransportDensity)
transportDomainDisc:add(basolateralTransport)

-- There is intentionally no urea_inter/apical coupling.  Interstitial urea
-- couples to the membrane only through the Basolateral subset.
transportDomainDisc:add(transportBoundary)

util.solver.defaults.approxSpace = transportSpace
transportSolver = util.solver.CreateSolver({
    type = "newton",
    convCheck = {
        type = "standard",
        iterations = 100,
        absolute = 1e-15,
        reduction = 1e-15,
        verbose = true
    },
    linSolver = "superlu"
})

timeDisc = ThetaTimeStep(transportDomainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(transportSolver)
timeIntegrator:set_time_step(dt)

ureaSolution = GridFunction(transportSpace)
Interpolate(0.0, ureaSolution, "urea_lumen")
Interpolate(0.0, ureaSolution, "urea_membrane")
Interpolate(0.0, ureaSolution, "urea_inter")

--------------------------------------------------------------------------------
-- Urea-only ParaView output
--------------------------------------------------------------------------------

-- Each VTK writer sees only the compartment on which its selected component
-- has exactly one nodal DoF.  This avoids the previous component-6/zero-DoF
-- failure and writes the complete spatial field of each compartment.
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

function writeCorrectedUreaVTK(solution, step, time)
    vtkLumen:print_subsets(
        vtkPrefix .. "_urea_lumen",
        solution,
        "Lumen",
        step,
        time
    )
    vtkMembrane:print_subsets(
        vtkPrefix .. "_urea_membrane",
        solution,
        "Membrane",
        step,
        time
    )
    vtkInter:print_subsets(
        vtkPrefix .. "_urea_inter",
        solution,
        "Inter",
        step,
        time
    )

    outputTimes[step] = time
    if step > lastOutputStep then
        lastOutputStep = step
    end
end

writeCorrectedUreaVTK(ureaSolution, 0, 0.0)

vtkObserver = LuaCallbackObserver()
function correctedVTKCallback(step, time, currentDt)
    writeCorrectedUreaVTK(
        vtkObserver:get_current_solution(),
        step,
        time
    )
    return 1
end
vtkObserver:set_callback("correctedVTKCallback")
timeIntegrator:attach_observer(vtkObserver)

transportSolver:init(AssembledOperator(transportDomainDisc))
transportSolver:prepare(ureaSolution)

local solveStartTime = GetClockS()
timeIntegrator:apply(
    ureaSolution,
    endTime,
    ureaSolution,
    0.0
)

--------------------------------------------------------------------------------
-- ParaView PVD time-series files
--------------------------------------------------------------------------------

local function writeCompartmentPVD(prefix)
    if ProcRank() ~= 0 then
        return
    end

    local pvdFile = assert(
        io.open(prefix .. ".pvd", "w"),
        "Cannot create PVD file: " .. prefix .. ".pvd"
    )
    local basename = string.match(prefix, "([^/]+)$")
    local extension = NumProcs() > 1 and "pvtu" or "vtu"

    pvdFile:write('<?xml version="1.0"?>\n')
    pvdFile:write(
        '<VTKFile type="Collection" version="0.1" ' ..
        'byte_order="LittleEndian">\n'
    )
    pvdFile:write("  <Collection>\n")

    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            pvdFile:write(string.format(
                '    <DataSet timestep="%.17g" group="" part="0" ' ..
                'file="%s_t%04d.%s"/>\n',
                outputTimes[step],
                basename,
                step,
                extension
            ))
        end
    end

    pvdFile:write("  </Collection>\n")
    pvdFile:write("</VTKFile>\n")
    pvdFile:close()
end

writeCompartmentPVD(vtkPrefix .. "_urea_lumen")
writeCompartmentPVD(vtkPrefix .. "_urea_membrane")
writeCompartmentPVD(vtkPrefix .. "_urea_inter")

--------------------------------------------------------------------------------
-- Runtime information
--------------------------------------------------------------------------------

local solveElapsedTime = GetClockS() - solveStartTime
local totalElapsedTime = GetClockS() - totalStartTime

print(string.format(
    "Transport runtime: %.6f mins",
    solveElapsedTime / 60.0
))
print(string.format(
    "Total runtime:     %.6f mins",
    totalElapsedTime / 60.0
))
print("ParaView array: Urea_normalized")
print(string.format(
    "Physical urea [mol/mm^3] = Urea_normalized * %.8g",
    concentrationScale
))
print("Urea ParaView output:")
print("  " .. vtkPrefix .. "_urea_lumen.pvd")
print("  " .. vtkPrefix .. "_urea_membrane.pvd")
print("  " .. vtkPrefix .. "_urea_inter.pvd")
