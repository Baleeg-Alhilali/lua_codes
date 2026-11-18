--------------------------------------------------------------------------------
-- One-way coupled nephron flow and transient urea transport
--
-- 1. Solve the validated steady Stokes--Darcy model on the current projected
--    generalized3 mesh.
-- 2. Freeze the converged Stokes and Darcy velocity fields.
-- 3. Solve transient urea advection--diffusion in Lumen, Membrane and Inter.
-- 4. Couple concentrations with passive Apical and Basolateral leak transport.
--
-- The lumen initially has concentration one.  Membrane and interstitium start
-- at zero.  Only the lumen Inlet is held at the prescribed concentration.
--------------------------------------------------------------------------------

local scriptPath = debug.getinfo(1, "S").source
if scriptPath:sub(1, 1) == "@" then scriptPath = scriptPath:sub(2) end
if scriptPath:sub(1, 1) ~= "/" then
    scriptPath = CurrentWorkingDirectory() .. "/" .. scriptPath
end
local MODEL_ROOT = assert(scriptPath:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. scriptPath)
ChangeDirectory(MODEL_ROOT)

-- This script leaves the solved domain, flow solution, pressure scaling and
-- flow equation objects available for the transient transport block below.
ug_load_script("code/flow/nephron_stokes_darcy_clean.lua")

--------------------------------------------------------------------------------
-- Transport parameters
--------------------------------------------------------------------------------

local endTime = util.GetParamNumber("-endTime", 100.0)
local dt = util.GetParamNumber("-dt", 1.0)
local outputEverySteps = util.GetParamNumber("-outputEverySteps", 1)
local diffusion = util.GetParamNumber("-diffusion", 1.67e-3) -- mm^2/s
local initialLumenConcentration = util.GetParamNumber("-initial", 1.0)
local inletConcentration = util.GetParamNumber("-inletConcentration", 1.0)
local apicalPermeability = util.GetParamNumber(
    "-apicalPermeability", 3.78e-4)       -- mm/s
local basolateralPermeability = util.GetParamNumber(
    "-basolateralPermeability", 3.70e-4)  -- mm/s
local membraneTransportDensity = util.GetParamNumber(
    "-membraneTransportDensity", 1.0)
local transportOutputPrefix = outputName .. "_urea"

assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0 and dt <= endTime, "-dt must be in (0, endTime]")
assert(outputEverySteps >= 1 and
       outputEverySteps == math.floor(outputEverySteps),
       "-outputEverySteps must be a positive integer")
assert(diffusion > 0.0, "-diffusion must be positive")
assert(initialLumenConcentration >= 0.0,
       "-initial must be non-negative")
assert(inletConcentration >= 0.0,
       "-inletConcentration must be non-negative")
assert(apicalPermeability >= 0.0,
       "-apicalPermeability must be non-negative")
assert(basolateralPermeability >= 0.0,
       "-basolateralPermeability must be non-negative")

print("=== TRANSIENT UREA ADVECTION--DIFFUSION AND MEMBRANE LEAK ===")
print(string.format("Time interval: 0 to %.17g s", endTime))
print(string.format("Time step: %.17g s", dt))
print(string.format("Initial lumen concentration: %.17g",
                    initialLumenConcentration))
print(string.format("Inlet concentration: %.17g", inletConcentration))

--------------------------------------------------------------------------------
-- Frozen flow coefficients and transport approximation space
--------------------------------------------------------------------------------

local solvedStokesVelocity = ExplicitGridFunctionVector(solution, "u,v,w")
local solvedScaledDarcyGradient = ExplicitGridFunctionGradient(
    solution, "p_darcy_scaled")

local transportMembraneVelocity = DarcyVelocityLinker()
transportMembraneVelocity:set_permeability(
    membranePermeability * pressureScale)
transportMembraneVelocity:set_viscosity(kinematicViscosity)
transportMembraneVelocity:set_density(1.0)
transportMembraneVelocity:set_gravity(ConstUserVector(0.0))
transportMembraneVelocity:set_pressure_gradient(solvedScaledDarcyGradient)

local transportInterVelocity = DarcyVelocityLinker()
transportInterVelocity:set_permeability(interPermeability * pressureScale)
transportInterVelocity:set_viscosity(kinematicViscosity)
transportInterVelocity:set_density(1.0)
transportInterVelocity:set_gravity(ConstUserVector(0.0))
transportInterVelocity:set_pressure_gradient(solvedScaledDarcyGradient)

transportSpace = ApproximationSpace(dom)
transportSpace:add_fct(
    "urea_lumen", "Lagrange", 1,
    "Lumen,Apical,Inlet,Outlet,TripleJunction")
transportSpace:add_fct(
    "urea_membrane", "Lagrange", 1,
    "Membrane,Apical,Basolateral,MembraneWall,TripleJunction")
-- TripleJunction does not touch an Inter volume on generalized3.  Excluding it
-- avoids the empty shadow rows that make a nonzero interface leak singular.
transportSpace:add_fct(
    "urea_inter", "Lagrange", 1,
    "Inter,Basolateral,OuterWall")
transportSpace:init_levels()
transportSpace:init_top_surface()
transportSpace:print_statistic("component")

--------------------------------------------------------------------------------
-- Transient advection--diffusion and passive interface transport
--------------------------------------------------------------------------------

transportDisc = DomainDiscretization(transportSpace)

local lumenTransport = ConvectionDiffusionFV1("urea_lumen", "Lumen")
lumenTransport:set_diffusion(diffusion)
lumenTransport:set_velocity(solvedStokesVelocity)
lumenTransport:set_upwind(FullUpwind())
transportDisc:add(lumenTransport)

local membraneTransport = ConvectionDiffusionFV1(
    "urea_membrane", "Membrane")
membraneTransport:set_diffusion(diffusion)
membraneTransport:set_velocity(transportMembraneVelocity)
membraneTransport:set_upwind(FullUpwind())
transportDisc:add(membraneTransport)

local interTransport = ConvectionDiffusionFV1("urea_inter", "Inter")
interTransport:set_diffusion(diffusion)
interTransport:set_velocity(transportInterVelocity)
interTransport:set_upwind(FullUpwind())
transportDisc:add(interTransport)

local apicalLeak = Leak({"urea_lumen", "urea_membrane"})
apicalLeak:set_permeability(apicalPermeability)
local apicalTransport = MembraneTransportFV1("Apical", apicalLeak)
apicalTransport:set_density_function(membraneTransportDensity)
transportDisc:add(apicalTransport)

local basolateralLeak = Leak({"urea_membrane", "urea_inter"})
basolateralLeak:set_permeability(basolateralPermeability)
local basolateralTransport = MembraneTransportFV1(
    "Basolateral", basolateralLeak)
basolateralTransport:set_density_function(membraneTransportDensity)
transportDisc:add(basolateralTransport)

local concentrationBoundary = DirichletBoundary()
concentrationBoundary:add(
    inletConcentration, "urea_lumen", "Inlet")
transportDisc:add(concentrationBoundary)
transportDisc:add(OneSideP1Constraints())

transportSolution = GridFunction(transportSpace)
Interpolate(initialLumenConcentration, transportSolution, "urea_lumen")
Interpolate(0.0, transportSolution, "urea_membrane")
Interpolate(0.0, transportSolution, "urea_inter")
transportDisc:adjust_solution(transportSolution)

local initialLumenMass = Integral(
    transportSolution, "urea_lumen", "Lumen")
local initialMembraneMass = Integral(
    transportSolution, "urea_membrane", "Membrane")
local initialInterMass = Integral(
    transportSolution, "urea_inter", "Inter")

--------------------------------------------------------------------------------
-- Scalable transient transport solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = transportSpace
local transportSolverConfig = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 6,
        absolute = 1.0e-14, reduction = 1.0e-10, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg", approxSpace = transportSpace, baseLevel = 0,
            smoother = {
                type = "gs", damping = 0.8, overlap = true,
                consistentInterfaces = true
            },
            cycle = "V", preSmooth = 3, postSmooth = 3,
            rap = true, transfer = "std", baseSolver = "superlu"
        },
        convCheck = {
            type = "standard", iterations = 200,
            absolute = 1.0e-14, reduction = 1.0e-10, verbose = true
        }
    }
}
local transportSolver = util.solver.CreateSolver(transportSolverConfig)

--------------------------------------------------------------------------------
-- Compartment time-series output: concentration and local velocity
--------------------------------------------------------------------------------

local outputTimes = {}
local lastOutputStep = 0
local compartmentOutputs = {
    {name = "lumen", subset = "Lumen", field = "urea_lumen",
     velocity = solvedStokesVelocity, writer = VTKOutput()},
    {name = "membrane", subset = "Membrane", field = "urea_membrane",
     velocity = transportMembraneVelocity, writer = VTKOutput()},
    {name = "inter", subset = "Inter", field = "urea_inter",
     velocity = transportInterVelocity, writer = VTKOutput()},
}
for _, compartment in ipairs(compartmentOutputs) do
    compartment.writer:clear_selection()
    compartment.writer:select(compartment.field, "Urea")
    compartment.writer:select_element(
        compartment.velocity, "velocity_mm_s")
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
    file:write('    <PCellData>\n')
    file:write('      <PDataArray type="Float32" Name="velocity_mm_s" ',
               'NumberOfComponents="3"/>\n')
    file:write('    </PCellData>\n')
    file:write(string.format('    <Piece Source="%s_t%04d.vtu"/>\n',
                             basename, step))
    file:write('  </PUnstructuredGrid>\n')
    file:write('</VTKFile>\n')
    file:close()
end

local function writeCompartmentVTK(currentSolution, step, time)
    for _, compartment in ipairs(compartmentOutputs) do
        local prefix = transportOutputPrefix .. "_" .. compartment.name
        compartment.writer:print_subsets(
            prefix, currentSolution, compartment.subset, step, time)
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

writeCompartmentVTK(transportSolution, 0, 0.0)

local timeDisc = ThetaTimeStep(transportDisc, 1.0)
local timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(transportSolver)
timeIntegrator:set_time_step(dt)
transportSolver:init(AssembledOperator(transportDisc))
transportSolver:prepare(transportSolution)

local vtkObserver = LuaCallbackObserver()
function coupledUreaVTKCallback(step, time, currentDt)
    if step % outputEverySteps == 0 or time >= endTime - 1.0e-12 then
        writeCompartmentVTK(vtkObserver:get_current_solution(), step, time)
    end
    return 1
end
vtkObserver:set_callback("coupledUreaVTKCallback")
timeIntegrator:attach_observer(vtkObserver)
timeIntegrator:apply(
    transportSolution, endTime, transportSolution, 0.0)

for _, compartment in ipairs(compartmentOutputs) do
    writeCompartmentPVD(
        transportOutputPrefix .. "_" .. compartment.name)
end

local finalLumenMass = Integral(
    transportSolution, "urea_lumen", "Lumen")
local finalMembraneMass = Integral(
    transportSolution, "urea_membrane", "Membrane")
local finalInterMass = Integral(
    transportSolution, "urea_inter", "Inter")
local initialTotalMass = initialLumenMass + initialMembraneMass + initialInterMass
local finalTotalMass = finalLumenMass + finalMembraneMass + finalInterMass

print(string.format("Initial masses: lumen=%.17e membrane=%.17e inter=%.17e",
    initialLumenMass, initialMembraneMass, initialInterMass))
print(string.format("Final masses: lumen=%.17e membrane=%.17e inter=%.17e",
    finalLumenMass, finalMembraneMass, finalInterMass))
print(string.format("Total mass change (inlet/outlet included): %.17e",
    finalTotalMass - initialTotalMass))
print("Transport output prefix: " .. transportOutputPrefix)
print("COUPLED_FLOW_TRANSIENT_UREA: COMPLETE")
