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
gridS = "third_trace_mm"
gridName = "runs/" .. gridS .. ".ugx"

numRefs = util.GetParamNumber("-numRefs", 1)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

endTime = util.GetParamNumber("-endTime", 2.0) -- [s]
dt = util.GetParamNumber("-dt", 0.1)           -- [s]

--------------------------------------------------------------------------------
-- Physical parameters
--------------------------------------------------------------------------------

-- The mesh coordinates are in micrometres, so diffusion is in um^2/s.
D_urea = 1.38
D_na = 1.33
D_cl = 2.03

-- Concentrations are in mM. Values imposed at the inlet can be overridden from
-- the command line, e.g. -ureaIn 1.0 -naIn 140.0 -clIn 140.0.
ureaIn = util.GetParamNumber("-ureaIn", 1.0)
naIn = util.GetParamNumber("-naIn", 140.0)
clIn = util.GetParamNumber("-clIn", 140.0)

-- Permeabilities are converted from m/s to um/s to match the mesh length unit.
m_to_um = 1.0e3
P_urea_apical = 1.0e-6 * m_to_um
P_urea_basolateral = 5.0e-7 * m_to_um
P_na_apical = 1.0e-7 * m_to_um
P_na_basolateral = 1.0e-7 * m_to_um
P_cl_apical = 1.0e-7 * m_to_um
P_cl_basolateral = 1.0e-7 * m_to_um

temperature = 310.15       -- [K]
voltageDrop = -70.0e-3     -- [V], defined as target potential - source potential
sourcePotential = 0.0       -- [V]
targetPotential = voltageDrop

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

-- Scalar algebra is required because subset-restricted functions give a
-- different number of unknowns at vertices on different subsets.
InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
     "OuterWall", "Inlet", "Outlet"},
    true)

balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

approxSpace = ApproximationSpace(dom)

local species = {"urea", "na", "cl"}
for _, s in ipairs(species) do
    approxSpace:add_fct(s .. "_lumen", "Lagrange", 1,
                        "Lumen,Apical,Inlet,Outlet")
    approxSpace:add_fct(s .. "_membrane", "Lagrange", 1,
                        "Membrane,Apical,Basolateral")
    approxSpace:add_fct(s .. "_inter", "Lagrange", 1,
                        "Inter,Basolateral,OuterWall")
end

approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "xy")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

dirichletBND = DirichletBoundary()
dirichletBND:add(ureaIn, "urea_lumen", "Inlet")
dirichletBND:add(naIn, "na_lumen", "Inlet")
dirichletBND:add(clIn, "cl_lumen", "Inlet")

domainDisc = DomainDiscretization(approxSpace)

local diffusion = {urea = D_urea, na = D_na, cl = D_cl}
for _, s in ipairs(species) do
    local lumenDisc = ConvectionDiffusionFV1(s .. "_lumen", "Lumen")
    lumenDisc:set_diffusion(diffusion[s])
    domainDisc:add(lumenDisc)

    local membraneDisc = ConvectionDiffusionFV1(s .. "_membrane", "Membrane")
    membraneDisc:set_diffusion(diffusion[s])
    domainDisc:add(membraneDisc)

    local interDisc = ConvectionDiffusionFV1(s .. "_inter", "Inter")
    interDisc:set_diffusion(diffusion[s])
    domainDisc:add(interDisc)
end

-- Neutral urea: ordinary concentration-driven passive transport.
ureaApicalLeak = Leak({"urea_lumen", "urea_membrane"})
ureaApicalLeak:set_permeability(P_urea_apical)
ureaApicalTransport = MembraneTransportFV1("Apical", ureaApicalLeak)
ureaApicalTransport:set_density_function(1.0)
domainDisc:add(ureaApicalTransport)

ureaBasolateralLeak = Leak({"urea_membrane", "urea_inter"})
ureaBasolateralLeak:set_permeability(P_urea_basolateral)
ureaBasolateralTransport = MembraneTransportFV1("Basolateral", ureaBasolateralLeak)
ureaBasolateralTransport:set_density_function(1.0)
domainDisc:add(ureaBasolateralTransport)

-- Charged solutes: four Leak inputs activate the GHK voltage-dependent flux.
-- Input order is {source concentration, target concentration,
--                 source potential, target potential}.
-- Empty potential names are replaced by fixed constants using zero-based input
-- indices 2 and 3. Thus phi_target - phi_source = -70 mV.
function addChargedTransport(interface, sourceFct, targetFct, permeability,
                             valency, objectPrefix)
    local ghkLeak = Leak({sourceFct, targetFct, "", ""})
    ghkLeak:set_constant(2, sourcePotential)
    ghkLeak:set_constant(3, targetPotential)
    ghkLeak:set_permeability(permeability)
    ghkLeak:set_temperature(temperature)
    ghkLeak:set_valency(valency)

    local transport = MembraneTransportFV1(interface, ghkLeak)
    transport:set_density_function(1.0)
    domainDisc:add(transport)

    -- Keep references alive and make them available for interactive inspection.
    _G[objectPrefix .. "Leak"] = ghkLeak
    _G[objectPrefix .. "Transport"] = transport
end

addChargedTransport("Apical", "na_lumen", "na_membrane",
                    P_na_apical, 1, "naApical")
addChargedTransport("Basolateral", "na_membrane", "na_inter",
                    P_na_basolateral, 1, "naBasolateral")
addChargedTransport("Apical", "cl_lumen", "cl_membrane",
                    P_cl_apical, -1, "clApical")
addChargedTransport("Basolateral", "cl_membrane", "cl_inter",
                    P_cl_basolateral, -1, "clBasolateral")

domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
-- MPI-friendly iterative solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace

solverDesc = {
    type = "newton",
    convCheck = {
        type = "standard", 
        iterations = 8, 
        absolute = 1e-12,
        reduction = 1e-12, 
        verbose = true
    },
    linSolver = {
        type = "bicgstab",
        precond = {
            type = "gmg",
            smoother = {
                type = "gs", overlap = true},
                cycle = "V", preSmooth = 4, postSmooth = 4,
                rap = true, baseLevel = numPreRefs,
            baseSolver = {
                type = "bicgstab",
                precond = {type = "ilu", overlap = true},
                convCheck = {
                    type = "standard", iterations = 200, absolute = 1e-10,
                    reduction = 0.5e-8, verbose = false
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
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)

u = GridFunction(approxSpace)
for _, s in ipairs(species) do
    Interpolate(0.0, u, s .. "_lumen")
    Interpolate(0.0, u, s .. "_membrane")
    Interpolate(0.0, u, s .. "_inter")
end

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(u, endTime, u, 0.0)
solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(string.format("Time integration runtime: %.6f mins", solve_elapsed_time / 60.0))
print(string.format("Total script runtime:      %.6f mins", total_elapsed_time / 60.0))
