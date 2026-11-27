--------------------------------------------------------------------------------
-- Refined nephron transient urea diffusion with passive membrane leak
--
-- The lumen starts at concentration one, while membrane and interstitium start
-- at zero. Only the lumen Inlet is held at one. Separate concentration fields
-- prevent artificial continuity across the membrane; passive Leak transport
-- couples Lumen -> Membrane at Apical and Membrane -> Inter at Basolateral.
--------------------------------------------------------------------------------

local scriptPath = debug.getinfo(1, "S").source
if scriptPath:sub(1, 1) == "@" then scriptPath = scriptPath:sub(2) end
if scriptPath:sub(1, 1) ~= "/" then
    scriptPath = CurrentWorkingDirectory() .. "/" .. scriptPath
end
local MODEL_ROOT = assert(scriptPath:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. scriptPath)
ChangeDirectory(MODEL_ROOT)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
AssertPluginsLoaded({"ConvectionDiffusion", "neuro_collection"})

dim = 3
local numRefs = util.GetParamNumber("-numRefs", 1)
assert(numRefs >= 1 and numRefs == math.floor(numRefs),
    "This validation requires an integer -numRefs >= 1")

local gridName = util.GetParam(
    "-grid", "mesh_generation/meshes/legacy_transport_meshes/third_trace_generalized3.ugx",
    "current coarse generalized-nephron NeuriteProjector mesh")
local endTime = util.GetParamNumber("-endTime", 0.1)
local dt = util.GetParamNumber("-dt", 0.1)
local outputEverySteps = util.GetParamNumber("-outputEverySteps", 1)
local diffusion = util.GetParamNumber("-diffusion", 1.67e-3) -- mm^2/s
local initialConcentration = util.GetParamNumber("-initial", 1.0)
local inletConcentration = util.GetParamNumber("-inlet", 1.0)
local apicalPermeability = util.GetParamNumber(
    "-apicalPermeability", 3.78e-4)       -- mm/s
local basolateralPermeability = util.GetParamNumber(
    "-basolateralPermeability", 3.70e-4)  -- mm/s
local membraneTransportDensity = util.GetParamNumber(
    "-membraneTransportDensity", 1.0)
local solverMode = util.GetParam("-solver", "gmg",
    "gmg for scalable runs; superlu only for small diagnostics")
local outputPrefix = util.GetParam(
    "-output",
    "Results_final/transport_diffusion_generalized3_r1_2026-11-12",
    "new VTK output prefix")

assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0 and dt <= endTime, "-dt must be in (0, endTime]")
assert(outputEverySteps >= 1 and outputEverySteps == math.floor(outputEverySteps),
    "-outputEverySteps must be a positive integer")
assert(diffusion > 0.0, "-diffusion must be positive")
assert(initialConcentration >= 0.0, "-initial must be non-negative")
assert(inletConcentration >= 0.0, "-inlet must be non-negative")
assert(apicalPermeability >= 0.0,
    "-apicalPermeability must be non-negative")
assert(basolateralPermeability >= 0.0,
    "-basolateralPermeability must be non-negative")
assert(membraneTransportDensity >= 0.0,
    "-membraneTransportDensity must be non-negative")
assert(solverMode == "gmg" or solverMode == "superlu",
    "-solver must be gmg or superlu")

local function fileContains(path, needle)
    local file = assert(io.open(path, "r"), "Cannot open grid file: " .. path)
    for line in file:lines() do
        if line:find(needle, 1, true) then
            file:close()
            return true
        end
    end
    file:close()
    return false
end

assert(fileContains(gridName, '<projection_handler'),
    "Mesh has no serialized ProjectionHandler: " .. gridName)
assert(fileContains(gridName, '<projector type="neurite"'),
    "Mesh has no serialized NeuriteProjector: " .. gridName)
assert(fileContains(gridName, '<vertex_attachment name="npSurfParams"'),
    "Mesh has no npSurfParams attachment: " .. gridName)

InitUG(dim, AlgebraType("CPU", 1))

-- The first load registers the custom projector attachment type.
local registrationDomain = Domain()
registrationDomain:create_additional_subset_handler("projSH")
LoadDomain(registrationDomain, gridName)
registrationDomain = nil
collectgarbage("collect")

dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {
    "Lumen", "Membrane", "Inter", "Apical", "Basolateral",
    "Inlet", "Outlet", "OuterWall", "TripleJunction", "MembraneWall"
}), "The mesh is missing a required current nephron subset")

-- Distribute the serialized coarse mesh before refinement, matching the clean
-- flow baseline. The current UG4 build transfers npSurfParams, so each rank's
-- NeuriteProjector refinement is valid and level 0 remains a usable GMG base.
balancer.parametersParsed = true
balancer.firstDistLvl = 0
balancer.RefineAndRebalanceDomain(dom, numRefs)

approxSpace = ApproximationSpace(dom)
local lumenFieldSubsets = "Lumen,Apical,Inlet,Outlet,TripleJunction"
local membraneFieldSubsets =
    "Membrane,Apical,Basolateral,MembraneWall,TripleJunction"
local interFieldSubsets =
    "Inter,Basolateral,OuterWall"
approxSpace:add_fct("urea_lumen", "Lagrange", 1, lumenFieldSubsets)
approxSpace:add_fct("urea_membrane", "Lagrange", 1, membraneFieldSubsets)
approxSpace:add_fct("urea_inter", "Lagrange", 1, interFieldSubsets)
approxSpace:init_levels()
approxSpace:init_top_surface()

print("=== REFINED TRANSIENT UREA LEAK TRANSPORT ===")
print("Grid: " .. gridName)
print("Refinements: " .. numRefs)
print("MPI processes: " .. NumProcs())
print(dom:domain_info():to_string())
approxSpace:print_statistic("component")

domainDisc = DomainDiscretization(approxSpace)
lumenDiffusion = ConvectionDiffusionFV1("urea_lumen", "Lumen")
lumenDiffusion:set_diffusion(diffusion)
membraneDiffusion = ConvectionDiffusionFV1("urea_membrane", "Membrane")
membraneDiffusion:set_diffusion(diffusion)
interDiffusion = ConvectionDiffusionFV1("urea_inter", "Inter")
interDiffusion:set_diffusion(diffusion)
domainDisc:add(lumenDiffusion)
domainDisc:add(membraneDiffusion)
domainDisc:add(interDiffusion)

apicalLeak = Leak({"urea_lumen", "urea_membrane"})
apicalLeak:set_permeability(apicalPermeability)
apicalTransport = MembraneTransportFV1("Apical", apicalLeak)
apicalTransport:set_density_function(membraneTransportDensity)
domainDisc:add(apicalTransport)

basolateralLeak = Leak({"urea_membrane", "urea_inter"})
basolateralLeak:set_permeability(basolateralPermeability)
basolateralTransport = MembraneTransportFV1(
    "Basolateral", basolateralLeak)
basolateralTransport:set_density_function(membraneTransportDensity)
domainDisc:add(basolateralTransport)

-- Apply the standard hanging-node constraint on the globally refined
-- projector hierarchy.  In particular, this keeps constrained interface
-- values from appearing as independent algebraic unknowns.
domainDisc:add(OneSideP1Constraints())

inletBoundary = DirichletBoundary()
inletBoundary:add(inletConcentration, "urea_lumen", "Inlet")
domainDisc:add(inletBoundary)
print(string.format("Initial lumen concentration: %.17g", initialConcentration))
print(string.format("Inlet Dirichlet concentration: %.17g", inletConcentration))
print(string.format("Apical leak permeability: %.17g mm/s", apicalPermeability))
print(string.format(
    "Basolateral leak permeability: %.17g mm/s", basolateralPermeability))

u = GridFunction(approxSpace)
Interpolate(initialConcentration, u, "urea_lumen")
Interpolate(0.0, u, "urea_membrane")
Interpolate(0.0, u, "urea_inter")
domainDisc:adjust_solution(u)

local initialLumenMass = Integral(u, "urea_lumen", "Lumen")
local initialMembraneMass = Integral(u, "urea_membrane", "Membrane")
local initialInterMass = Integral(u, "urea_inter", "Inter")
local lumenVolume = Integral(1.0, u, "Lumen")
print(string.format("Initial lumen mass: %.17e", initialLumenMass))
print(string.format("Lumen volume: %.17e mm^3", lumenVolume))
print(string.format("Initial membrane mass: %.17e", initialMembraneMass))
print(string.format("Initial inter mass: %.17e", initialInterMass))

util.solver.defaults.approxSpace = approxSpace
local solverConfig = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 6,
        absolute = 1.0e-14, reduction = 1.0e-10, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg", approxSpace = approxSpace, baseLevel = 0,
            smoother = {
                type = "gs", damping = 0.8, overlap = true,
                consistentInterfaces = true
            },
            cycle = "V", preSmooth = 3, postSmooth = 3,
            rap = true, transfer = "std",
            -- The coupled three-compartment operator is poorly conditioned on
            -- the small level-0 grid.  A direct solve is confined to that GMG
            -- coarse level; the refined system still uses BiCGStab + GMG.
            baseSolver = "superlu"
        },
        convCheck = {
            type = "standard", iterations = 150,
            absolute = 1.0e-14, reduction = 1.0e-10, verbose = true
        }
    }
}
if solverMode == "superlu" then
    solverConfig.linSolver = "superlu"
end
print("Solver mode: " .. solverMode)
solver = util.solver.CreateSolver(solverConfig)

-- ParaView output is compartment-specific and time-aware.  Each compartment
-- gets one PVD collection for the full transient.  Every PVD entry points to
-- a PVTU master for one time step.  UG4 writes a VTU directly when NumProcs()
-- is one, so create the corresponding one-piece PVTU master explicitly for
-- mpirun -np 1 runs.
local outputTimes = {}
local lastOutputStep = 0
local compartmentOutputs = {
    {name = "lumen", subset = "Lumen", field = "urea_lumen",
     writer = VTKOutput()},
    {name = "membrane", subset = "Membrane", field = "urea_membrane",
     writer = VTKOutput()},
    {name = "inter", subset = "Inter", field = "urea_inter",
     writer = VTKOutput()},
}
for _, compartment in ipairs(compartmentOutputs) do
    compartment.writer:select(compartment.field, "Urea")
end

local function writeSinglePiecePVTU(prefix, step, time)
    if ProcRank() ~= 0 or NumProcs() ~= 1 then return end

    local fileName = string.format("%s_t%04d.pvtu", prefix, step)
    local file = assert(io.open(fileName, "w"),
        "Cannot create PVTU file: " .. fileName)
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n')
    file:write('<VTKFile type="PUnstructuredGrid" version="0.1" ',
               'byte_order="LittleEndian">\n')
    file:write(string.format('  <Time timestep="%.17g"/>\n', time))
    file:write('  <PUnstructuredGrid GhostLevel="0">\n')
    file:write('    <PPoints>\n')
    file:write('      <PDataArray type="Float32" ',
               'NumberOfComponents="3"/>\n')
    file:write('    </PPoints>\n')
    file:write('    <PPointData>\n')
    file:write('      <PDataArray type="Float32" Name="Urea" ',
               'NumberOfComponents="1"/>\n')
    file:write('    </PPointData>\n')
    file:write(string.format('    <Piece Source="%s_t%04d.vtu"/>\n',
                             basename, step))
    file:write('  </PUnstructuredGrid>\n')
    file:write('</VTKFile>\n')
    file:close()
end

local function writeCompartmentVTK(solution, step, time)
    for _, compartment in ipairs(compartmentOutputs) do
        local prefix = outputPrefix .. "_" .. compartment.name
        compartment.writer:print_subsets(
            prefix, solution, compartment.subset, step, time)
        writeSinglePiecePVTU(prefix, step, time)
    end
    outputTimes[step] = time
    if step > lastOutputStep then lastOutputStep = step end
end

local function writeCompartmentPVD(prefix)
    if ProcRank() ~= 0 then return end

    local file = assert(io.open(prefix .. ".pvd", "w"),
        "Cannot create PVD file: " .. prefix .. ".pvd")
    local basename = assert(prefix:match("([^/]+)$"))
    file:write('<?xml version="1.0"?>\n')
    file:write('<VTKFile type="Collection" version="0.1" ',
               'byte_order="LittleEndian">\n')
    file:write('  <Collection>\n')
    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            file:write(string.format(
                '    <DataSet timestep="%.17g" group="" part="0" ' ..
                'file="%s_t%04d.pvtu"/>\n',
                outputTimes[step], basename, step))
        end
    end
    file:write('  </Collection>\n')
    file:write('</VTKFile>\n')
    file:close()
end

writeCompartmentVTK(u, 0, 0.0)
timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
solver:init(AssembledOperator(domainDisc))
solver:prepare(u)
vtkObserver = LuaCallbackObserver()
function diffusionVTKCallback(step, time, currentDt)
    if step % outputEverySteps == 0 or time >= endTime - 1.0e-12 then
        writeCompartmentVTK(
            vtkObserver:get_current_solution(), step, time)
    end
    return 1
end
vtkObserver:set_callback("diffusionVTKCallback")
timeIntegrator:attach_observer(vtkObserver)
timeIntegrator:apply(u, endTime, u, 0.0)
for _, compartment in ipairs(compartmentOutputs) do
    writeCompartmentPVD(outputPrefix .. "_" .. compartment.name)
end

local finalLumenMass = Integral(u, "urea_lumen", "Lumen")
local finalMembraneMass = Integral(u, "urea_membrane", "Membrane")
local finalInterMass = Integral(u, "urea_inter", "Inter")
local initialTotalMass =
    initialLumenMass + initialMembraneMass + initialInterMass
local finalTotalMass = finalLumenMass + finalMembraneMass + finalInterMass

print(string.format("Final lumen mass: %.17e", finalLumenMass))
print(string.format("Final membrane mass: %.17e", finalMembraneMass))
print(string.format("Final inter mass: %.17e", finalInterMass))
print(string.format("Lumen mass change: %.17e",
    finalLumenMass - initialLumenMass))
print(string.format("Membrane mass change: %.17e",
    finalMembraneMass - initialMembraneMass))
print(string.format("Inter mass change: %.17e",
    finalInterMass - initialInterMass))
print(string.format("Total mass change (includes inlet supply): %.17e",
    finalTotalMass - initialTotalMass))
print("VTK output prefix: " .. outputPrefix)
print("TRANSIENT_UREA_APICAL_BASOLATERAL_LEAK: COMPLETE")
print("Bounds and compartment exchange require independent VTK inspection.")
