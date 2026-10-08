-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--------------------------------------------------------------------------------
-- COMPLETE 2D UREA / Na+ / Cl- / K+ MEMBRANE-TRANSPORT MODEL
--
-- Run:
--   ~/UG4/ug4/bin/ugshell -ex code/transport/2D_trans_combined_print.lua
--
-- Mechanisms:
--   Urea: passive Leak (neutral).
--   Na+:  voltage-dependent GHK Leak, valency +1.
--   Cl-:  voltage-dependent GHK Leak, valency -1.
--   K+:   voltage-dependent GHK Leak, valency +1.
--   Na/K pump at Basolateral: 3 Na+ Membrane -> Inter and
--                              2 K+ Inter -> Membrane per ATP cycle.
--
-- Command-line parameters and defaults:
--   -numRefs 4                    mesh refinements
--   -numPreRefs 0                 refinements before load balancing
--   -endTime 15.0                 final time [s]
--   -dt 0.1                       time step [s]
--   -diffusionUrea 0.138          urea diffusion [um^2/s]
--   -diffusionNa 0.133            Na+ diffusion [um^2/s]
--   -diffusionCl 0.203            Cl- diffusion [um^2/s]
--   -diffusionK 0.196             K+ diffusion [um^2/s]
--   -inletConcentration 1.0       common LeftBnd value [mM]
--   -apicalPermeability 1e-6      common apical permeability [m/s]
--   -basolateralPermeability 5e-7 common basolateral permeability [m/s]
--   -pumpMaxFlux 1e-22            maximum cycles [mol/(pump s)]
--   -pumpKNa 10.0                 Na half-saturation [mM]
--   -pumpKK 1.5                   K half-saturation [mM]
--   -pumpDensity 500              pumps/um^2
--   -initialNaMembrane 10.0       initial cell/membrane Na [mM]
--
-- Example using a deliberately larger pump rate for visualization:
--   ~/UG4/ug4/bin/ugshell -ex code/transport/2D_trans_combined_print.lua \
--       -pumpMaxFlux 1e-19 -pumpDensity 500
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 2
gridS = "2D_trans2"
gridName = "ProMeshFiles/runs/" .. gridS .. ".ugx"

numRefs = util.GetParamNumber("-numRefs", 4)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
vtk_file_name = "results/" .. gridS .. "_combined_" .. numRefs
endTime = util.GetParamNumber("-endTime", 15.0) -- [s]
dt = util.GetParamNumber("-dt", 0.1)           -- [s]

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- Each species has its own diffusion coefficient. The mesh coordinates are in
-- micrometres, so all four values must be supplied in um^2/s.
D_urea = util.GetParamNumber("-diffusionUrea", 1.38e-1)
D_na = util.GetParamNumber("-diffusionNa", 1.33e-1)
D_cl = util.GetParamNumber("-diffusionCl", 2.03e-1)
D_k = util.GetParamNumber("-diffusionK", 1.96e-1)

-- Lookup table used by the species loop below. The table keys must match the
-- strings in the species list: "urea", "na", "cl", and "k".
diffusionCoefficient = {
    urea = D_urea,
    na = D_na,
    cl = D_cl,
    k = D_k
}

-- All inlet concentrations are identical and expressed in mM.
concentration_left = util.GetParamNumber("-inletConcentration", 1.0)

-- All species have the same permeability at a given interface. Permeabilities
-- use um/s to match the mesh length unit.
m_to_um = 1.0e6
P_apical = util.GetParamNumber("-apicalPermeability", 1.0e-6) * m_to_um
P_basolateral = util.GetParamNumber("-basolateralPermeability", 5.0e-7) * m_to_um

temperature = 310.15   -- [K]

-- Fixed electrical potential drop in the positive transport direction.
-- The charged Leak implementation uses v = phi_target - phi_source.
phi_source = 0.0
phi_target = -70.0e-3  -- [V], so phi_target - phi_source = -70 mV

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

-- Subset-restricted functions require scalar block algebra.
InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {
        "Lumen", "Membrane", "Basolateral", "Apical", "Inter",
        "Toplum", "Topmem", "Topint", "LeftBnd", "RightBnd"
    },
    true
)

-- util.CreateDomain(gridFile, preRefinements, requiredSubsets, distribute)
--   gridFile:        UGX geometry path.
--   preRefinements:  refinements performed before parallel distribution.
--   requiredSubsets: names that must exist in the UGX subset handler.
--   distribute:      true enables parallel domain distribution.

balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

-- Each species has a separate function in each compartment. This permits a
-- concentration jump across Apical and Basolateral instead of enforcing
-- artificial continuity through the membrane.
approxSpace = ApproximationSpace(dom)

species = {"urea", "na", "cl", "k"}
for _, s in ipairs(species) do
    -- add_fct(name, basis, order, subsets)
    -- creates a scalar first-order Lagrange concentration field only on the
    -- listed subsets. Separate compartment fields allow concentration jumps.
    approxSpace:add_fct(
        s .. "_lumen", "Lagrange", 1,
        "Lumen,Apical,LeftBnd,Toplum"
    )
    approxSpace:add_fct(
        s .. "_membrane", "Lagrange", 1,
        "Membrane,Apical,Basolateral,Topmem"
    )
    approxSpace:add_fct(
        s .. "_inter", "Lagrange", 1,
        "Inter,Basolateral,RightBnd,Topint"
    )
end

approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "x")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

dirichletBND = DirichletBoundary()
-- DirichletBoundary:add(value, functionName, subsetName) fixes one field to a
-- constant value on a boundary. Here LeftBnd supplies solute and RightBnd is
-- the interstitial reservoir; extracellular K is fixed to 1 mM.
dirichletBND:add(concentration_left, "urea_lumen", "LeftBnd")
dirichletBND:add(concentration_left, "na_lumen", "LeftBnd")
dirichletBND:add(concentration_left, "cl_lumen", "LeftBnd")
dirichletBND:add(concentration_left, "k_lumen", "LeftBnd")
dirichletBND:add(0.0, "urea_inter", "RightBnd")
dirichletBND:add(0.0, "na_inter", "RightBnd")
dirichletBND:add(0.0, "cl_inter", "RightBnd")
dirichletBND:add(1.0, "k_inter", "RightBnd")

domainDisc = DomainDiscretization(approxSpace)

for _, s in ipairs(species) do
    -- ConvectionDiffusionFV1(functionName, volumeSubset) constructs the FV1
    -- volume equation. diffusionCoefficient[s] selects the coefficient that
    -- belongs to the current species s.
    lumenDiffusion = ConvectionDiffusionFV1(s .. "_lumen", "Lumen")
    lumenDiffusion:set_diffusion(diffusionCoefficient[s])
    domainDisc:add(lumenDiffusion)

    membraneDiffusion = ConvectionDiffusionFV1(
        s .. "_membrane", "Membrane"
    )
    membraneDiffusion:set_diffusion(diffusionCoefficient[s])
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
ureaApicalLeak:set_permeability(P_apical)
ureaApicalTransport = MembraneTransportFV1("Apical", ureaApicalLeak)
ureaApicalTransport:set_density_function(0.4)
domainDisc:add(ureaApicalTransport)

ureaBasolateralLeak = Leak({"urea_membrane", "urea_inter"})
ureaBasolateralLeak:set_permeability(P_basolateral)
ureaBasolateralTransport = MembraneTransportFV1("Basolateral", ureaBasolateralLeak)
ureaBasolateralTransport:set_density_function(1.0)
domainDisc:add(ureaBasolateralTransport)

--------------------------------------------------------------------------------
-- Voltage-dependent Na+ and Cl- membrane transport
--------------------------------------------------------------------------------

-- Four-input Leak evaluates the Goldman-Hodgkin-Katz flux. In this installed
-- Neuro Collection version the voltage mode is selected in the constructor
-- only when at least one potential slot contains a function name. Therefore an
-- existing interface function is supplied as a harmless activation anchor;
-- set_constant then overrides both potential inputs with the fixed voltages.
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
    transporter:set_permeability(permeability)
    transporter:set_temperature(temperature)
    transporter:set_valency(valency)

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
    P_apical, 1, 0.4, "naApical"
)
addChargedTransport(
    "Basolateral", "na_membrane", "na_inter", "urea_membrane",
    P_basolateral, 1, 1.0, "naBasolateral"
)

-- Cl- has valency -1, reversing its electrical response relative to Na+.
addChargedTransport(
    "Apical", "cl_lumen", "cl_membrane", "urea_lumen",
    P_apical, -1, 0.4, "clApical"
)
addChargedTransport(
    "Basolateral", "cl_membrane", "cl_inter", "urea_membrane",
    P_basolateral, -1, 1.0, "clBasolateral"
)

-- K+ passive electrodiffusion through the same GHK formulation as Na+.
addChargedTransport(
    "Apical", "k_lumen", "k_membrane", "urea_lumen",
    P_apical, 1, 0.4, "kApical"
)
addChargedTransport(
    "Basolateral", "k_membrane", "k_inter", "urea_membrane",
    P_basolateral, 1, 1.0, "kBasolateral"
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
            cycle = "V", preSmooth = 4, postSmooth = 4,
            rap = true, baseLevel = numPreRefs,
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

u = GridFunction(approxSpace)
for _, s in ipairs(species) do
    Interpolate(0.0, u, s .. "_lumen")
    Interpolate(0.0, u, s .. "_membrane")
    Interpolate(0.0, u, s .. "_inter")
end
-- The standard pump requires intracellular Na and extracellular K. Without
-- these nonzero initial substrates its rate is exactly zero at startup.
Interpolate(util.GetParamNumber("-initialNaMembrane", 10.0),
            u, "na_membrane")
Interpolate(1.0, u, "k_inter")

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

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
    vtkMembrane:print_subsets(
        vtk_file_name .. "_membrane", uOut, "Membrane", step, time
    )
    vtkInter:print_subsets(
        vtk_file_name .. "_inter", uOut, "Inter", step, time
    )
end

writeVTK(u, 0, 0.0)

vtkObserver = LuaCallbackObserver()

function vtkCallback(step, time, currentDt)
    writeVTK(vtkObserver:get_current_solution(), step, time)
    return 1
end

vtkObserver:set_callback("vtkCallback")

--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
-- ThetaTimeStep(discretization, theta): theta=1 is implicit Euler.
timeIntegrator = SimpleTimeIntegrator(timeDisc)
-- apply(solutionOut, endTime, initialSolution, startTime) advances the model.
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
timeIntegrator:attach_observer(vtkObserver)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)

vtkLumen:write_time_pvd(vtk_file_name .. "_lumen", u)
vtkMembrane:write_time_pvd(vtk_file_name .. "_membrane", u)
vtkInter:write_time_pvd(vtk_file_name .. "_inter", u)

solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time / 60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time / 60.0))
