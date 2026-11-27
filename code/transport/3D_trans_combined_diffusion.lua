-- Resolve all relative paths from the Model directory.
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

local totalStartTime = GetClockS()

--------------------------------------------------------------------------------
-- Pure-diffusion 3D Urea / Na+ / Cl- / K+ membrane-transport model
--
-- There is no prescribed velocity and no flow equation. Each solute diffuses
-- in Lumen, Membrane, and Inter. Separate compartment fields permit jumps at
-- Apical and Basolateral, where the same membrane mechanisms as the 2D model
-- are assembled.
--------------------------------------------------------------------------------

dim = 3
local refinementLevel = util.GetParamNumber("-numRefs", 1)
assert(refinementLevel >= 0 and refinementLevel == math.floor(refinementLevel),
    "-numRefs must be a non-negative integer")

-- Always start from the original coarse mesh. Its serialized
-- NeuriteProjectors generate the curved geometry on every new level.
local gridName = util.GetParam(
    "-grid", "mesh_generation/meshes/legacy_transport_meshes/third_trace_mm_transport.ugx",
    "coarse transport projector UGX mesh")
local numPreRefs = 0

local gridFile = io.open(gridName, "r")
assert(gridFile ~= nil, "Cannot open grid file: " .. gridName)
gridFile:close()

local endTime = util.GetParamNumber("-endTime", 2.0)
local dt = util.GetParamNumber("-dt", 0.1)
assert(endTime > 0.0, "-endTime must be positive")
assert(dt > 0.0, "-dt must be positive")

local enableVTK = util.GetParamBool("-vtk", true)
local timeTag = tostring(endTime):gsub("%.", "p")
local vtkPrefix = util.GetParam(
    "-vtkPrefix",
    string.format("Results_final/third_trace_3d_diffusion_T%s_ref%d",
                  timeTag, refinementLevel),
    "VTK/PVTU output prefix"
)

--------------------------------------------------------------------------------
-- Concentrations and transport parameters (mesh length unit: um)
--------------------------------------------------------------------------------

-- Independent inlet values make physiological concentrations possible while
-- retaining the 2D model's common 1 mM default.
local inlet = {
    urea = util.GetParamNumber("-ureaIn", 1.0),
    na = util.GetParamNumber("-naIn", 1.0),
    cl = util.GetParamNumber("-clIn", 1.0),
    k = util.GetParamNumber("-kIn", 1.0)
}

-- Free-solution diffusion coefficients [um^2/s].
local diffusion = {
    urea = util.GetParamNumber("-diffusionUrea", 1.80e3),
    na = util.GetParamNumber("-diffusionNa", 1.33e3),
    cl = util.GetParamNumber("-diffusionCl", 2.60e3),
    k = util.GetParamNumber("-diffusionK", 2.50e3)
}
local membraneHindrance = util.GetParamNumber("-membraneHindrance", 1.0)

-- Match the 2D interface permeabilities after converting m/s to um/s.
local mToUm = 1.0e6
local pApical = util.GetParamNumber("-apicalPermeability", 1.0e-6) * mToUm
local pBasolateral = util.GetParamNumber(
    "-basolateralPermeability", 5.0e-7) * mToUm
local apicalDensity = util.GetParamNumber("-apicalDensity", 0.4)
local basolateralDensity = util.GetParamNumber("-basolateralDensity", 1.0)
local enableUreaTransport = util.GetParamBool("-ureaTransport", true)
local enableChargedTransport = util.GetParamBool("-chargedTransport", true)
local enablePump = util.GetParamBool("-pump", true)
local diffusionSpeciesCount = util.GetParamNumber("-diffusionSpeciesCount", 4)
assert(diffusionSpeciesCount >= 1 and diffusionSpeciesCount <= 4 and
       diffusionSpeciesCount == math.floor(diffusionSpeciesCount),
       "-diffusionSpeciesCount must be an integer from 1 through 4")

local temperature = 310.15
local phiSource = 0.0
local phiTarget = -70.0e-3

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())
print("Grid: " .. gridName)
print("VTK output prefix: " .. vtkPrefix)

-- Register the serialized neurite attachment type before the actual load.
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
}), "The mesh is missing a required 3D transport subset")

-- Projector refinement is created on one rank before the completed hierarchy
-- is distributed. This preserves an identical curved hierarchy in MPI runs.
balancer.firstDistLvl = refinementLevel
balancer.RefineAndRebalanceDomain(dom, refinementLevel - numPreRefs)
-- Project while the coordinates still match the serialized millimetre
-- centerline, then scale every generated level to micrometres. This retains
-- valid curved refinement and gives the transport operator natural um units.
ScaleDomain(dom, 1.0e3, 1.0e3, 1.0e3)

print("Domain info:")
print(dom:domain_info():to_string())

local species = {"urea", "na", "cl", "k"}
approxSpace = ApproximationSpace(dom)
for _, s in ipairs(species) do
    approxSpace:add_fct(
        s .. "_lumen", "Lagrange", 1,
        "Lumen,Apical,Inlet,Outlet,TripleJunction"
    )
    approxSpace:add_fct(
        s .. "_membrane", "Lagrange", 1,
        "Membrane,Apical,Basolateral,TripleJunction"
    )
    approxSpace:add_fct(
        s .. "_inter", "Lagrange", 1,
        "Inter,Basolateral,OuterWall,Inlet,Outlet,TripleJunction"
    )
end
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()
-- Keep UG4's native ordering. OrderLex cannot consistently order these
-- restricted compartment fields and damages the level-smoother layout.

--------------------------------------------------------------------------------
-- Diffusion and boundary conditions
--------------------------------------------------------------------------------

domainDisc = DomainDiscretization(approxSpace)

-- The inlet is the only prescribed concentration boundary. Outlet and
-- OuterWall retain the natural zero-diffusive-flux condition.
local inletBoundary = DirichletBoundary()
for speciesIndex, s in ipairs(species) do
    inletBoundary:add(inlet[s], s .. "_lumen", "Inlet")

    if speciesIndex <= diffusionSpeciesCount then
        local lumenDisc = ConvectionDiffusionFV1(s .. "_lumen", "Lumen")
        lumenDisc:set_diffusion(diffusion[s])
        domainDisc:add(lumenDisc)

        local membraneDisc = ConvectionDiffusionFV1(s .. "_membrane", "Membrane")
        membraneDisc:set_diffusion(diffusion[s] * membraneHindrance)
        domainDisc:add(membraneDisc)

        local interDisc = ConvectionDiffusionFV1(s .. "_inter", "Inter")
        interDisc:set_diffusion(diffusion[s])
        domainDisc:add(interDisc)
    end
end

--------------------------------------------------------------------------------
-- Apical and basolateral membrane fluxes (same mechanisms as the 2D model)
--------------------------------------------------------------------------------

if enableUreaTransport then
    local ureaApicalLeak = Leak({"urea_lumen", "urea_membrane"})
    ureaApicalLeak:set_permeability(pApical)
    local ureaApicalTransport = MembraneTransportFV1("Apical", ureaApicalLeak)
    ureaApicalTransport:set_density_function(apicalDensity)
    domainDisc:add(ureaApicalTransport)

    local ureaBasolateralLeak = Leak({"urea_membrane", "urea_inter"})
    ureaBasolateralLeak:set_permeability(pBasolateral)
    local ureaBasolateralTransport = MembraneTransportFV1(
        "Basolateral", ureaBasolateralLeak)
    ureaBasolateralTransport:set_density_function(basolateralDensity)
    domainDisc:add(ureaBasolateralTransport)
end

local function addChargedTransport(
    interface, sourceFct, targetFct, voltageAnchorFct,
    permeability, valency, density, objectPrefix)

    local leak = Leak({sourceFct, targetFct, voltageAnchorFct, ""})
    leak:set_constant(2, phiSource)
    leak:set_constant(3, phiTarget)
    leak:set_permeability(permeability)
    leak:set_temperature(temperature)
    leak:set_valency(valency)

    local disc = MembraneTransportFV1(interface, leak)
    disc:set_density_function(density)
    domainDisc:add(disc)

    _G[objectPrefix .. "Leak"] = leak
    _G[objectPrefix .. "Transport"] = disc
end

if enableChargedTransport then
    addChargedTransport("Apical", "na_lumen", "na_membrane", "urea_lumen",
                        pApical, 1, apicalDensity, "naApical")
    addChargedTransport("Basolateral", "na_membrane", "na_inter",
                        "urea_membrane", pBasolateral, 1,
                        basolateralDensity, "naBasolateral")
    addChargedTransport("Apical", "cl_lumen", "cl_membrane", "urea_lumen",
                        pApical, -1, apicalDensity, "clApical")
    addChargedTransport("Basolateral", "cl_membrane", "cl_inter",
                        "urea_membrane", pBasolateral, -1,
                        basolateralDensity, "clBasolateral")
    addChargedTransport("Apical", "k_lumen", "k_membrane", "urea_lumen",
                        pApical, 1, apicalDensity, "kApical")
    addChargedTransport("Basolateral", "k_membrane", "k_inter",
                        "urea_membrane", pBasolateral, 1,
                        basolateralDensity, "kBasolateral")
end

-- Basolateral Na/K ATPase: 3 Na_membrane -> Na_inter and
-- 2 K_inter -> K_membrane per cycle.
if enablePump then
    local naKPump = NaKPump({
        "na_membrane", "na_inter", "k_membrane", "k_inter"
    })
    naKPump:set_max_flux(util.GetParamNumber("-pumpMaxFlux", 1.0e-22))
    naKPump:set_k_na(util.GetParamNumber("-pumpKNa", 10.0))
    naKPump:set_k_k(util.GetParamNumber("-pumpKK", 1.5))
    -- mol/(um^2 s) -> mM*um/s.
    naKPump:set_scale_fluxes({1.0e15, 1.0e15})
    local naKPumpTransport = MembraneTransportFV1("Basolateral", naKPump)
    naKPumpTransport:set_density_function(
        util.GetParamNumber("-pumpDensity", 500.0))
    domainDisc:add(naKPumpTransport)
end

domainDisc:add(inletBoundary)

--------------------------------------------------------------------------------
-- Initial conditions
--------------------------------------------------------------------------------

u = GridFunction(approxSpace)
for _, s in ipairs(species) do
    -- The entire lumen starts at its fixed inlet concentration.
    Interpolate(inlet[s], u, s .. "_lumen")
    Interpolate(0.0, u, s .. "_membrane")
    Interpolate(0.0, u, s .. "_inter")
end

-- Preserve the 2D pump's optional initial substrates. Set either to zero on
-- the command line if a completely solute-free non-lumen initial state is
-- desired.
Interpolate(util.GetParamNumber("-initialNaMembrane", 10.0),
            u, "na_membrane")
Interpolate(util.GetParamNumber("-initialKInter", 1.0),
            u, "k_inter")

--------------------------------------------------------------------------------
-- Solver and time integration
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace
local solverDesc = {
    type = "newton",
    convCheck = {
        type = "standard", iterations = 8, absolute = 1e-14,
        reduction = 1e-8, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg",
            approxSpace = approxSpace,
            baseLevel = 0,
            -- Coupled ILU is more robust here than scalar Gauss--Seidel.
            smoother = {
                type = "ilu",
                damping = 0.7,
                overlap = true,
                consistentInterfaces = true
            },
            cycle = "V",
            preSmooth = 4,
            postSmooth = 4,
            rap = false,
            transfer = "std",
            -- Keep even the level-0 solve iterative and MPI-compatible. This
            -- avoids the memory growth and serial bottleneck of SuperLU.
            baseSolver = {
                type = "bicgstab",
                precond = {
                    type = "ilu",
                    damping = 0.7,
                    inversionEps = 1e-8,
                    overlap = true,
                    consistentInterfaces = true
                },
                convCheck = {
                    type = "standard", iterations = 500,
                    absolute = 1e-14, reduction = 1e-9,
                    verbose = false
                }
            }
        },
        convCheck = {
            type = "standard", iterations = 200,
            absolute = 1e-14, reduction = 1e-9,
            verbose = true
        }
    }
}
solver = util.solver.CreateSolver(solverDesc)

--------------------------------------------------------------------------------
-- Combined ParaView output
--------------------------------------------------------------------------------

local outputTimes = {}
local lastOutputStep = 0
local vtkLumen, vtkMembrane, vtkInter, vtkObserver

local function writeCombinedPVTU(step, time)
    if ProcRank() ~= 0 then return end
    local fileName = string.format("%s_t%04d.pvtu", vtkPrefix, step)
    local file = assert(io.open(fileName, "w"),
        "Cannot create combined PVTU: " .. fileName)
    local baseName = assert(vtkPrefix:match("([^/]+)$"))
    local compartments = {"lumen", "membrane", "inter"}

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
    for _, label in ipairs({"Urea", "Na", "Cl", "K"}) do
        file:write(string.format(
            '      <PDataArray type="Float32" Name="%s" ' ..
            'NumberOfComponents="1"/>\n', label))
    end
    file:write('    </PPointData>\n')
    for _, compartment in ipairs(compartments) do
        if NumProcs() == 1 then
            file:write(string.format(
                '    <Piece Source="%s_%s_t%04d.vtu"/>\n',
                baseName, compartment, step))
        else
            for rank = 0, NumProcs() - 1 do
                file:write(string.format(
                    '    <Piece Source="%s_%s_p%d_t%04d.vtu"/>\n',
                    baseName, compartment, rank, step))
            end
        end
    end
    file:write('  </PUnstructuredGrid>\n')
    file:write('</VTKFile>\n')
    file:close()

    outputTimes[step] = time
    if step > lastOutputStep then lastOutputStep = step end
end

local function writeVTK(solution, step, time)
    vtkLumen:print_subsets(vtkPrefix .. "_lumen", solution,
                           "Lumen", step, time)
    vtkMembrane:print_subsets(vtkPrefix .. "_membrane", solution,
                              "Membrane", step, time)
    vtkInter:print_subsets(vtkPrefix .. "_inter", solution,
                           "Inter", step, time)
    writeCombinedPVTU(step, time)
end

if enableVTK then
    vtkLumen = VTKOutput()
    vtkMembrane = VTKOutput()
    vtkInter = VTKOutput()
    for _, s in ipairs(species) do
        local label = s == "na" and "Na"
            or (s == "cl" and "Cl" or (s == "k" and "K" or "Urea"))
        vtkLumen:select(s .. "_lumen", label)
        vtkMembrane:select(s .. "_membrane", label)
        vtkInter:select(s .. "_inter", label)
    end
    writeVTK(u, 0, 0.0)

    vtkObserver = LuaCallbackObserver()
    function transport3DVTKCallback(step, time, currentDt)
        writeVTK(vtkObserver:get_current_solution(), step, time)
        return 1
    end
    vtkObserver:set_callback("transport3DVTKCallback")
end

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
if enableVTK then timeIntegrator:attach_observer(vtkObserver) end

local solveStartTime = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)

print(string.format("Time integration runtime: %.6f mins",
                    (GetClockS() - solveStartTime) / 60.0))
print(string.format("Total script runtime:      %.6f mins",
                    (GetClockS() - totalStartTime) / 60.0))
