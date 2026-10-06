-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 3
-- This transport-safe copy gives shared junction edges dedicated subsets,
-- removes detached OuterWall quadrilaterals, and gives the lumen-only inlet
-- and outlet face interiors their own subsets. These classifications prevent
-- refinement from creating missing or orphan compartment DoFs.
gridS = "third_trace_ogrid"
gridName = "runs/"..gridS..".ugx"

-- The base UGX is deliberately coarse and has a block-like cross-section.
-- At least one refinement is needed for the embedded NeuriteProjectors to add
-- projected surface vertices and reveal the cylindrical lumen. Override with
-- -numRefs N when a finer (N > 1) or deliberately coarse (N = 0) grid is wanted.
numRefs = util.GetParamNumber("-numRefs", 0)
numPreRefs = 0 -- 0

-- Command-line overrides make short output checks possible, for example:
--   ugshell -ex code/Dif_itterative_multi.lua -endTime 0.1 -dt 0.1
endTime = util.GetParamNumber("-endTime", 2.0) -- [s]
dt = util.GetParamNumber("-dt", 0.1)          -- [s]
-- One switch controls the complete output pipeline. Run with "-vtk false" to
-- disable VTK creation without leaving an initial write or PVD call active.
enableVTK = util.GetParamBool("-vtk", true)

-- A literal decimal point in the base name is interpreted by UG4 as the start
-- of a file extension. Replace it so runs such as endTime=0.1 do not collapse
-- to the same truncated name and overwrite the three compartment outputs.
endTimeTag = string.gsub(tostring(endTime), "%.", "p")
vtk_file_name = "Results/" .. gridS .. "_T" .. endTimeTag .. "_ref" .. numRefs
print("VTK output prefix: " .. vtk_file_name)

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- The mesh coordinates are in micrometres, so diffusion is in um^2/s.
D_na = 1.33e3
D_k = 2.5e3
D_cl = 2.6e3
D_urea = 1.8e3
-- Hindrance factor and effective diffusion inside the cell layer.
hindrance_factor = 1.0e-3
D_na_cell = D_na * hindrance_factor
D_k_cell = D_k * hindrance_factor
D_cl_cell = D_cl * hindrance_factor
D_urea_cell = D_urea * hindrance_factor

diffusionCoefficient = {
    urea = D_urea,
    na = D_na,
    cl = D_cl,
    k = D_k
}
--------------------------------------------------------------------------------
-- Apical membrane permeabilities [m/s]
--------------------------------------------------------------------------------

P_na_apical = 1.0e-6
P_k_apical = 5.0e-6
P_cl_apical = 8.0e-7
P_urea_apical = 1.0e-6

-- GHK parameters. Leak uses voltage = target potential - source potential.
-- These values therefore impose a -70 mV drop in the positive transport
-- direction for every charged membrane transporter.
temperature = 310.15       -- [K]
phi_source = 0.0           -- [V]
phi_target = -70.0e-3      -- [V]
--------------------------------------------------------------------------------
-- Basolateral membrane permeabilities [m/s]
--------------------------------------------------------------------------------

P_na_basolateral = 5.0e-7
P_k_basolateral = 2.0e-5
P_cl_basolateral = 1.0e-6
P_urea_basolateral = 1.0e-6

--------------------------------------------------------------------------------
-- Na/K pump parameters
--------------------------------------------------------------------------------

Imax_pump = 10.0           -- [A/m^2], converted from 1 mA/cm^2
K_na_pump = 30.0           -- [mol/m^3], converted from 0.03 mmol/cm^3
K_k_pump = 27.0            -- [mol/m^3], converted from 0.027 mmol/cm^3

--------------------------------------------------------------------------------
-- Inlet concentrations [mol/m^3]
--------------------------------------------------------------------------------

C_in_na = 140.0
C_in_k = 4.5
C_in_cl = 110.0
C_in_urea = 5.0

--------------------------------------------------------------------------------
-- Initial lumen concentrations [mol/m^3]
--------------------------------------------------------------------------------

C_initial_na_lumen = 140.0
C_initial_k_lumen = 4.5
C_initial_cl_lumen = 110.0
C_initial_urea_lumen = 5.0

--------------------------------------------------------------------------------
-- Initial cell concentrations [mol/m^3]
--------------------------------------------------------------------------------

C_initial_na_cell = 15.0
C_initial_k_cell = 130.0
C_initial_cl_cell = 20.0
C_initial_urea_cell = 5.0



--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

-- The UGX stores its NeuriteProjectors in a serialized ProjectionHandler.
-- Its auxiliary subset handler ("projSH") must exist before LoadDomain is
-- called; util.CreateDomain does not create it and would therefore leave the
-- model with ordinary linear refinement. Linear refinement preserves the
-- coarse square cross-section instead of projecting new Apical/Basolateral
-- vertices onto the cylindrical neurite geometry.
dom = Domain()
dom:create_additional_subset_handler("projSH")
print("Loading the embedded Apical and Basolateral neurite projectors ...")
LoadDomain(dom, gridName)

-- neuro_collection supplies this compatibility loader. It installs the
-- serialized UGX NeuriteProjector as the domain refinement projector without
-- requiring any modification to UG4 core. It must be called after LoadDomain
-- (so npSurfParams exists) and before the first refinement.
InstallNeuriteProjectorFromUGX(dom, gridName)

assert(
    util.CheckSubsets(
        dom,
        {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
         "OuterWall", "Inlet", "Outlet", "LumenInletFace",
         "LumenOutletFace"}
    ),
    "The domain is missing one or more required subsets"
)

-- If serial pre-refinement is requested later, it must happen only after the
-- projector has been loaded. Every new boundary vertex will then be projected.
if numPreRefs > 0 then
    local preRefiner = GlobalDomainRefiner(dom)
    for i = 1, numPreRefs do
        preRefiner:refine()
    end
    delete(preRefiner)
end

-- Refine and distribute the grid among all MPI processes. The loaded UGX
-- ProjectionHandler is used automatically; do not replace it with a new
-- ProjectionHandler, since that would discard the two NeuriteProjectors.
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)
species = {"urea", "na", "cl", "k"} 
for _, s in ipairs(species) do
    -- add_fct(name, basis, order, subsets)
    -- creates a scalar first-order Lagrange concentration field only on the
    -- listed subsets. Separate compartment fields allow concentration jumps.
    approxSpace:add_fct(
        s .. "_lumen", "Lagrange", 1,
        "Lumen,Apical,Inlet,Outlet,LumenInletFace,LumenOutletFace," ..
        "TripleJunction"
    )
    approxSpace:add_fct(
        s .. "_membrane", "Lagrange", 1,
        -- Membrane hexahedra also touch the end-cap vertex subsets.  FV1
        -- requires one DoF at every corner of every Membrane hexahedron.
        "Membrane,Apical,Basolateral,Inlet,Outlet," ..
        "BasolateralJunction,TripleJunction"
    )
    approxSpace:add_fct(
        s .. "_inter", "Lagrange", 1,
        "Inter,Basolateral,OuterWall,Inlet,Outlet," ..
        "BasolateralJunction,TripleJunction"
    )
end
approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

-- Order the DoFs:
OrderLex (approxSpace, "xy")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

dirichletBND = DirichletBoundary()
--inlet concentration 
dirichletBND:add(C_in_na, "na_lumen", "Inlet")
dirichletBND:add(C_in_cl, "cl_lumen", "Inlet")
dirichletBND:add(C_in_k, "k_lumen", "Inlet")
dirichletBND:add(C_in_urea, "urea_lumen", "Inlet")
dirichletBND:add(C_in_na, "na_lumen", "LumenInletFace")
dirichletBND:add(C_in_cl, "cl_lumen", "LumenInletFace")
dirichletBND:add(C_in_k, "k_lumen", "LumenInletFace")
dirichletBND:add(C_in_urea, "urea_lumen", "LumenInletFace")
--outlet boundry condetion 
dirichletBND:add(0.0, "na_lumen", "Outlet")
dirichletBND:add(0.0, "cl_lumen", "Outlet")
dirichletBND:add(0.0, "k_lumen", "Outlet")
dirichletBND:add(0.0, "urea_lumen", "Outlet")
dirichletBND:add(0.0, "na_lumen", "LumenOutletFace")
dirichletBND:add(0.0, "cl_lumen", "LumenOutletFace")
dirichletBND:add(0.0, "k_lumen", "LumenOutletFace")
dirichletBND:add(0.0, "urea_lumen", "LumenOutletFace")
--outer wll boundry condetion 
dirichletBND:add(0.0, "na_inter", "OuterWall")
dirichletBND:add(0.0, "cl_inter", "OuterWall")
dirichletBND:add(C_in_k, "k_inter", "OuterWall")
dirichletBND:add(0.0, "urea_inter", "OuterWall")
--------------------------------------------------------------------------------
-- Domain Eqs
--------------------------------------------------------------------------------

-- Create the container before adding any volume equations, membrane transport
-- equations, or boundary conditions to it.
domainDisc = DomainDiscretization(approxSpace)

for _, s in ipairs(species) do
    -- ConvectionDiffusionFV1(functionName, volumeSubset) constructs the FV1
    -- volume equation. set_diffusion(D) supplies the scalar diffusion value.
    
    lumenDiffusion = ConvectionDiffusionFV1(s .. "_lumen", "Lumen")
    lumenDiffusion:set_diffusion(diffusionCoefficient[s])
    domainDisc:add(lumenDiffusion)

    membraneDiffusion = ConvectionDiffusionFV1(s .. "_membrane", "Membrane")
    membraneDiffusion:set_diffusion(diffusionCoefficient[s]*hindrance_factor)
    domainDisc:add(membraneDiffusion)

    interDiffusion = ConvectionDiffusionFV1(s .. "_inter", "Inter")
    interDiffusion:set_diffusion(diffusionCoefficient[s])
    domainDisc:add(interDiffusion)
end
--------------------------------------------------------------------------------
-- Neutral urea membrane transport
--------------------------------------------------------------------------------
-- Two-input Leak gives ordinary passive flux P * (c_source - c_target).
-- Leak({sourceConcentration, targetConcentration}) defines positive flux from
-- the first field to the second. set_permeability(P) supplies P in um/s here.
ureaApicalLeak = Leak({"urea_lumen", "urea_membrane"})
ureaApicalLeak:set_permeability(P_urea_apical)
ureaApicalTransport = MembraneTransportFV1("Apical", ureaApicalLeak)
ureaApicalTransport:set_density_function(0.4)
domainDisc:add(ureaApicalTransport)

ureaBasolateralLeak = Leak({"urea_membrane", "urea_inter"})
ureaBasolateralLeak:set_permeability(P_urea_basolateral)
ureaBasolateralTransport = MembraneTransportFV1("Basolateral", ureaBasolateralLeak)
ureaBasolateralTransport:set_density_function(1.0)
domainDisc:add(ureaBasolateralTransport)

--------------------------------------------------------------------------------
-- Voltage-dependent Na+ and Cl- membrane transport
--------------------------------------------------------------------------------
function addChargedTransport(
    interface, sourceFct, targetFct, voltageAnchorFct,
    permeability, valency, density, name
)
    -- Parameters:
    --   interface        interface subset on which flux is assembled
    --   sourceFct        source-side concentration field
    --   targetFct        target-side concentration field
    --   voltageAnchorFct existing field used to activate this UG4 GHK branch
    --   permeability     membrane permeability in um/s
    --   valency          signed ionic charge (+1 for Na/K, -1 for Cl)
    --   density          dimensionless transporter-area multiplier
    --   name             prefix used to retain globally inspectable objects
    --
    -- Leak's four inputs are {c_source, c_target, phi_source, phi_target}.
    transporter = Leak({sourceFct, targetFct, voltageAnchorFct, ""})
    -- set_constant(index, value) uses zero-based constructor-input indices.
    -- Indices 2 and 3 replace the two potential inputs with values in volts.
    transporter:set_constant(2, phi_source)
    transporter:set_constant(3, phi_target)
    transporter:set_temperature(temperature)
    transporter:set_valency(valency)
    transporter:set_permeability(permeability)

    transportDisc = MembraneTransportFV1(interface, transporter)
    -- MembraneTransportFV1(interfaceSubset, mechanism) converts the mechanism's
    -- flux into equal and opposite FV1 contributions on both interface sides.
    -- set_density_function(value) multiplies that flux by a constant density.
    transportDisc:set_density_function(density)
    domainDisc:add(transportDisc)

    -- Preserve Lua references and useful names for interactive inspection.
    _G[name .. "Leak"] = transporter
    _G[name .. "Transport"] = transportDisc
end

-- Na+ has valency +1.
addChargedTransport(
    "Apical", "na_lumen", "na_membrane", "urea_lumen",
    P_na_apical, 1, 0.4, "naApical"
)
addChargedTransport(
    "Basolateral", "na_membrane", "na_inter", "urea_membrane",
    P_na_basolateral, 1, 1.0, "naBasolateral"
)

-- Cl- has valency -1, reversing its electrical response relative to Na+.
addChargedTransport(
    "Apical", "cl_lumen", "cl_membrane", "urea_lumen",
    P_cl_apical, -1, 0.4, "clApical"
)
addChargedTransport(
    "Basolateral", "cl_membrane", "cl_inter", "urea_membrane",
    P_cl_basolateral, -1, 1.0, "clBasolateral"
)

-- K+ passive electrodiffusion through the same GHK formulation as Na+.
addChargedTransport(
    "Apical", "k_lumen", "k_membrane", "urea_lumen",
    P_k_apical, 1, 0.4, "kApical"
)
addChargedTransport(
    "Basolateral", "k_membrane", "k_inter", "urea_membrane",
    P_k_basolateral, 1, 1.0, "kBasolateral"
)

-- Active basolateral 3 Na+ out / 2 K+ in transport, superimposed on GHK.
-- NaKPump({Na_i, Na_o, K_i, K_o}) requires this exact field order.
-- Its kinetic rate is:
-- Jcycle = Jmax * Na_i^3/(KNa^3 + Na_i^3)
--               * K_o^2/(KK^2 + K_o^2).
naKPump = NaKPump({"na_membrane", "na_inter", "k_membrane", "k_inter"})
naKPump:set_max_flux(util.GetParamNumber("-pumpMaxFlux", 1.0e-22))
-- set_max_flux: maximum ATPase cycle flux [mol/(pump s)].
naKPump:set_k_na(util.GetParamNumber("-pumpKNa", 10.0))
-- set_k_na: intracellular Na half-saturation parameter [mM].
naKPump:set_k_k(util.GetParamNumber("-pumpKK", 1.5))
-- set_k_k: extracellular K half-saturation parameter [mM].
naKPump:set_scale_fluxes({1.0e15, 1.0e15})
-- One scale for each output (Na and K); converts mol/(um^2 s) into the
-- concentration-flux unit mM*um/s used by this micrometre geometry.
naKPumpTransport = MembraneTransportFV1("Basolateral", naKPump)
naKPumpTransport:set_density_function(util.GetParamNumber("-pumpDensity", 500.0))
domainDisc:add(naKPumpTransport)
domainDisc:add(dirichletBND)


--------------------------------------------------------------------------------
-- MPI-friendly iterative solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace

-- Geometric multigrid is needed for this large refined diffusion system.
-- The overlapping Gauss-Seidel smoother and iterative overlapping-ILU base
-- solver are MPI-compatible; no rank-local direct factorization is used.
solverDesc = {
    type = "newton",
    convCheck = {
        type = "standard", 
        iterations = 8, absolute = 1e-8,
        reduction = 1e-8, verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg",
            smoother = {type = "gs", overlap = true},
            cycle = "V",
            preSmooth = 4,
            postSmooth = 4,
            rap = true,
            baseLevel = numPreRefs,
            baseSolver = {
                type = "bicgstab",
                precond = {type = "ilu", overlap = true},
                convCheck = {
                    type = "standard", iterations = 200,
                    absolute = 1e-10, reduction = 0.5e-8, verbose = false
                }
            }
        },
        convCheck = {
            type = "standard", iterations = 100, absolute = 1e-8,
            reduction = 0.5e-8, verbose = true
        }
    }
}

solver = util.solver.CreateSolver(solverDesc)

--------------------------------------------------------------------------------
-- initial condetions
--------------------------------------------------------------------------------

u = GridFunction(approxSpace)
for _, s in ipairs(species) do
    Interpolate(0.0, u, s .. "_lumen")
    Interpolate(0.0, u, s .. "_membrane")
    Interpolate(0.0, u, s .. "_inter")
end
-- The standard pump requires intracellular Na and extracellular K. Without
-- these nonzero initial substrates its rate is exactly zero at startup.
Interpolate(10.0,u, "na_membrane")
Interpolate(C_in_k, u, "k_inter")
Interpolate(C_in_na,u, "na_lumen")
Interpolate(C_in_k,u, "k_lumen")

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

if enableVTK then
vtkLumen = VTKOutput()
vtkMembrane = VTKOutput()
vtkInter = VTKOutput()
for _, s in ipairs(species) do
    local label = s == "na" and "Na" or (s == "k" and "K" or (s == "cl" and "Cl" or "Urea"))
    vtkLumen:select(s .. "_lumen", label)
    vtkMembrane:select(s .. "_membrane", label)
    vtkInter:select(s .. "_inter", label)
end

function writeVTK(uOut, step, time)
    vtkLumen:print_subsets(
        vtk_file_name .. "_lumen", uOut, "Lumen", step, time
    )
    -- The Membrane subset contains volume hexahedra.  Its fields now include
    -- every vertex subset touched by those hexahedra, so nodal output is valid.
    vtkMembrane:print_subsets(
        vtk_file_name .. "_membrane", uOut, "Membrane", step, time
    )
    vtkInter:print_subsets(
        vtk_file_name .. "_inter", uOut, "Inter", step, time
    )
    outputTimes[step] = time
    if step > lastOutputStep then lastOutputStep = step end
end

-- Keep the exact observer times so that the PVD collection can refer to the
-- compartment files produced by print_subsets(). UG4's built-in
-- write_time_pvd() assumes that all restricted functions were printed as
-- eight separate geometry-subset files; that assumption is false here.
outputTimes = {}
lastOutputStep = 0
writeVTK(u, 0, 0.0)

vtkObserver = LuaCallbackObserver()

function vtkCallback(step, time, currentDt)
    writeVTK(vtkObserver:get_current_solution(), step, time)
    return 1
end

vtkObserver:set_callback("vtkCallback")
end -- complete VTK setup, including the initial write

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
-- ThetaTimeStep(discretization, theta): theta=1 is implicit Euler.
timeIntegrator = SimpleTimeIntegrator(timeDisc)
-- apply(solutionOut, endTime, initialSolution, startTime) advances the model.
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
if enableVTK then
    timeIntegrator:attach_observer(vtkObserver)
end

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)

if enableVTK then
function writeCompartmentPVD(prefix)
    -- Only one rank writes the shared collection index. With multiple MPI
    -- ranks, print_subsets() creates a PVTU master for each time; otherwise it
    -- creates one VTU file directly.
    if ProcRank() ~= 0 then return end

    local collection = assert(io.open(prefix .. ".pvd", "w"))
    local basename = string.match(prefix, "([^/]+)$")
    local extension = NumProcs() > 1 and "pvtu" or "vtu"

    collection:write('<?xml version="1.0"?>\n')
    collection:write('<VTKFile type="Collection" version="0.1">\n')
    collection:write('  <Collection>\n')
    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            collection:write(string.format(
                '    <DataSet timestep="%.17g" part="0" file="%s_t%04d.%s"/>\n',
                outputTimes[step], basename, step, extension
            ))
        end
    end
    collection:write('  </Collection>\n')
    collection:write('</VTKFile>\n')
    collection:close()
end

writeCompartmentPVD(vtk_file_name .. "_lumen")
writeCompartmentPVD(vtk_file_name .. "_membrane")
writeCompartmentPVD(vtk_file_name .. "_inter")
end -- final VTK time-series indices

solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time / 60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time / 60.0))
