-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- Physical constants
--------------------------------------------------------------------------------

Faraday = 96485.0          -- [C/mol]
gas_constant = 8.3144      -- [J/(mol K)] (R_i in the COMSOL model)
Temperature = 310.15       -- [K], body temperature
beta = Faraday / (gas_constant * Temperature) -- [1/V]
eps_ghk = 1.0e-12          -- Stabilizer used in GHK expressions

--------------------------------------------------------------------------------
-- Ionic valences
--------------------------------------------------------------------------------

zNa = 1.0
zK = 1.0
zCl = -1.0
zUrea = 0.0

--------------------------------------------------------------------------------
-- Diffusion coefficients [m^2/s]
--------------------------------------------------------------------------------

D_na = 1.7e-9
D_k = 2.5e-9
D_cl = 2.6e-9
D_urea = 1.8e-9

-- Hindrance factor and effective diffusion inside the cell layer.
hindrance_factor = 1.0e-3
D_na_cell = D_na * hindrance_factor
D_k_cell = D_k * hindrance_factor
D_cl_cell = D_cl * hindrance_factor
D_urea_cell = D_urea * hindrance_factor

--------------------------------------------------------------------------------
-- Apical membrane permeabilities [m/s]
--------------------------------------------------------------------------------

P_na_apical = 1.0e-6
P_k_apical = 5.0e-6
P_cl_apical = 8.0e-7
P_urea_apical = 1.0e-6

--------------------------------------------------------------------------------
-- Basolateral membrane permeabilities [m/s]
--------------------------------------------------------------------------------

P_na_basolateral = 5.0e-7
P_k_basolateral = 2.0e-5
P_cl_basolateral = 1.0e-6
P_urea_basolateral = 1.0e-6

--------------------------------------------------------------------------------
-- Paracellular permeabilities [m/s]
--------------------------------------------------------------------------------

P_na_paracellular = 3.0e-5
P_k_paracellular = 2.5e-5
P_cl_paracellular = 2.0e-5
P_urea_paracellular = 1.0e-6
paracellular_factor = 0.05

--------------------------------------------------------------------------------
-- Na/K pump parameters
--------------------------------------------------------------------------------

Imax_pump = 10.0           -- [A/m^2], converted from 1 mA/cm^2
K_na_pump = 30.0           -- [mol/m^3], converted from 0.03 mmol/cm^3
K_k_pump = 27.0            -- [mol/m^3], converted from 0.027 mmol/cm^3

--------------------------------------------------------------------------------
-- Geometry [m]
--------------------------------------------------------------------------------

tube_length = 1.0e-3
tube_diameter = 30.0e-6
cell_thickness = 20.0e-6
extracellular_thickness = 200.0e-6

--------------------------------------------------------------------------------
-- Inlet flow and bolus parameters
--------------------------------------------------------------------------------

inflow_ultrafiltrate = 10.0e-12 / 60.0 -- [m^3/s], converted from 10 nl/min
bolus_volume = 100.0e-12                -- [m^3], converted from 100 nl
initial_time = 0.0                      -- [s]

-- The final COMSOL history overwrites Q_in with an electrical expression.
-- Use the physical inlet flow for the UG4 model unless that overwrite is verified.
Q_in = inflow_ultrafiltrate             -- [m^3/s]
bolus_duration = bolus_volume / Q_in    -- [s]
bolus_end_time = initial_time + bolus_duration -- [s]

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

fixed_charge_concentration = 125.0 -- [mol/m^3], A_fixed in COMSOL
fixed_protein_concentration = 80.0 -- [mol/m^3], C_fixed_pro in COMSOL

--------------------------------------------------------------------------------
-- Reflection coefficients
--------------------------------------------------------------------------------

sigma_na = 0.75
sigma_k = 0.60
sigma_cl = 0.30
sigma_urea = 0.10
sigma_protein = 1.0

--------------------------------------------------------------------------------
-- Hydraulic permeability [m/(Pa s)]
--------------------------------------------------------------------------------

Lp_membrane = 1.0e-12
Lp_paracellular = 5.0e-8

--------------------------------------------------------------------------------
-- Hydrostatic pressures [Pa]
--------------------------------------------------------------------------------

mmHg_to_Pa = 133.322387415
pressure_lumen = 12.0 * mmHg_to_Pa
pressure_cell = 7.0 * mmHg_to_Pa
pressure_extracellular = 5.0 * mmHg_to_Pa

--------------------------------------------------------------------------------
-- Fluid properties
--------------------------------------------------------------------------------

density_blood = 1050.0             -- [kg/m^3]
density_ultrafiltrate = 1000.0     -- [kg/m^3]
density_water = 1000.0             -- [kg/m^3]
viscosity_blood = 0.7e-3           -- [Pa s], converted from 0.7 cP
viscosity_ultrafiltrate = 0.7e-3   -- [Pa s], converted from 0.7 cP

--------------------------------------------------------------------------------
-- Initial electrical coefficient retained from COMSOL
--------------------------------------------------------------------------------

-- COMSOL name: phi_lum_int. Despite its name, this expression has the form of
-- an ionic conductivity coefficient rather than an initial electric potential.
phi_lumen_initial_coefficient =
    Faraday^2 / (gas_constant * Temperature) *
    (zNa * D_na * C_initial_na_lumen +
     zK * D_k * C_initial_k_lumen +
     C_initial_cl_lumen * D_cl)

--------------------------------------------------------------------------------
-- UG4: steady Navier-Stokes flow in the nephron lumen
--------------------------------------------------------------------------------

PluginRequired("NavierStokes")

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")
ug_load_script("Examples/navier_stokes_util.lua")
ug_load_script("util/conv_rates_static.lua")

-- The ProMesh coordinates are currently used exactly as stored in the UGX
-- file. No domain scaling is applied at this stage.
dim = 3
gridName = util.GetParam(
    "-grid",
    "ProMeshFiles/fourth_trace_all_element_counts.ugx",
    "nephron grid")
numRefs = util.GetParamNumber("-numRefs", 0, "number of grid refinements")
numPreRefs = util.GetParamNumber(
    "-numPreRefs", 0, "refinements before MPI distribution")
-- FV1 is the default because UG4's explicit no-normal-stress outlet element
-- is available for FV discretizations (not for the standard FE variant).
discType = util.GetParam("-discType", "fv1", "Navier-Stokes discretization")
vorder = util.GetParamNumber("-vorder", 1, "velocity order")
porder = util.GetParamNumber("-porder", 1, "pressure order")

bStokes = util.HasParamOption("-stokes", "solve Stokes instead of Navier-Stokes")
bNoLaplace = util.HasParamOption("-nolaplace", "use deformation tensor formulation")
bExactJac = util.HasParamOption("-exactjac", "use the exact Jacobian")
bPecletBlend = util.HasParamOption("-pecletblend", "use Peclet blending")
upwind = util.GetParam("-upwind", "full", "upwind type")
stab = util.GetParam("-stab", "flow", "stabilization type")
diffLength = util.GetParam("-difflength", "cor", "diffusion length type")

lumen_subset = "Lumen"
inlet_subset = "Inlet"
wall_subset = "Apical"
outlet_subset = "Outlet"
inactive_subsets = "Membrane,Basolateral,OuterWall,Inter"

-- Retain the inlet speed obtained previously from 10 nl/min and the nominal
-- 30 micrometre lumen diameter.
lumen_radius = 0.5 * tube_diameter
lumen_inlet_area = math.pi * lumen_radius^2
inlet_speed = inflow_ultrafiltrate / lumen_inlet_area -- [m/s]

-- Inward unit vector calculated from the Inlet face and its adjacent Lumen
-- volume in fourth_trace_all_element_counts.ugx.
inlet_direction_x = -0.6711542847787513
inlet_direction_y =  0.7389250208227904
inlet_direction_z =  0.05951251654199889

function inletVelocity(x, y, z, t)
    return inlet_speed * inlet_direction_x,
           inlet_speed * inlet_direction_y,
           inlet_speed * inlet_direction_z
end

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

function CreateDomain()
    InitUG(dim, AlgebraType("CPU", 1))
    local requiredSubsets = {
        lumen_subset, inlet_subset, wall_subset, outlet_subset
    }
    return util.CreateAndDistributeDomain(
        gridName, numRefs, numPreRefs, requiredSubsets)
end

function CreateApproxSpace(dom)
    local approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder)
    approxSpace:init_top_surface()
    approxSpace:print_statistic()
    return approxSpace
end

--------------------------------------------------------------------------------
-- Stationary incompressible Navier-Stokes discretization
--------------------------------------------------------------------------------

-- The UGX coordinates are stored in micrometres and are intentionally not
-- rescaled. Convert the SI kinematic viscosity to the coefficient required
-- when one numerical coordinate unit represents one micrometre:
--     nu_mesh = nu_SI / (metres per mesh unit)
mesh_length_unit = 1.0e-6 -- [m/mesh unit]
kinematic_viscosity_si =
    viscosity_ultrafiltrate / density_ultrafiltrate -- [m^2/s]
kinematic_viscosity =
    kinematic_viscosity_si / mesh_length_unit -- [mesh unit m/s]

globalNSDisc = nil

function CreateDomainDisc(approxSpace)
    local fctCmp = approxSpace:names()
    local nsDisc = NavierStokes(fctCmp, {lumen_subset}, discType)

    nsDisc:set_exact_jacobian(bExactJac)
    nsDisc:set_stokes(bStokes)
    nsDisc:set_laplace(not bNoLaplace)
    nsDisc:set_kinematic_viscosity(kinematic_viscosity)
    globalNSDisc = nsDisc

    local actualPOrder = approxSpace:lfeid(dim):order()
    local actualVOrder = approxSpace:lfeid(0):order()

    if discType == "fv1" or discType == "fvcr" then
        nsDisc:set_upwind(upwind)
        nsDisc:set_peclet_blend(bPecletBlend)
    end

    if discType == "fv1" then
        nsDisc:set_stabilization(stab, diffLength)
        nsDisc:set_pac_upwind(true)
    end

    if discType == "fe" and actualPOrder == actualVOrder then
        nsDisc:set_stabilization(3)
    end

    if discType == "fe" or discType == "fv" then
        nsDisc:set_quad_order(math.pow(actualVOrder, dim) + 2)
    end

    -- Prescribed uniform velocity at the lumen inlet.
    local inletDisc = NavierStokesInflow(nsDisc)
    inletDisc:add("inletVelocity", inlet_subset)

    -- No-slip velocity condition on the lumen/apical interface.
    local wallDisc = NavierStokesWall(nsDisc)
    wallDisc:add(wall_subset)

    -- Explicit zero normal-stress boundary condition at the outlet.
    local outletDisc = NavierStokesNoNormalStressOutflow(nsDisc)
    outletDisc:add(outlet_subset)

    -- The UGX contains several non-fluid materials. The approximation-space
    -- helper creates degrees of freedom on the complete mesh, so explicitly
    -- constrain those otherwise equation-free unknowns. Without this block,
    -- the global Navier--Stokes system is singular.
    local inactiveDisc = DirichletBoundary()
    inactiveDisc:add(0.0, "u", inactive_subsets)
    inactiveDisc:add(0.0, "v", inactive_subsets)
    inactiveDisc:add(0.0, "w", inactive_subsets)
    inactiveDisc:add(0.0, "p", inactive_subsets)

    local domainDisc = DomainDiscretization(approxSpace)
    domainDisc:add(nsDisc)
    domainDisc:add(inletDisc)
    domainDisc:add(wallDisc)
    domainDisc:add(outletDisc)
    domainDisc:add(inactiveDisc)

    return domainDisc
end

--------------------------------------------------------------------------------
-- Nonlinear and linear solvers
--------------------------------------------------------------------------------

function CreateSolver(approxSpace)
    local base = LU()
    local smoother

    if discType == "fvcr" or discType == "fecr" then
        smoother = ComponentGaussSeidel(0.1, {"p"}, {1, 2}, {1})
    elseif discType == "fv1" then
        smoother = ILU()
        smoother:set_damp(0.7)
    elseif discType == "fe" and porder == vorder then
        smoother = ILU()
        smoother:set_damp(0.7)
    else
        smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
    end

    -- These command-line parsers and factories follow the UG4 example and
    -- allow solver/multigrid settings to be changed without editing this file.
    local smooth = util.smooth.parseParams()
    smoother = util.smooth.create(smooth)

    local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
    local gmg = util.gmg.create(
        approxSpace, smoother, numPreSmooth, numPostSmooth,
        cycle, base, baseLev, bRAP)
    gmg:add_prolongation_post_process(AverageComponent("p"))

    local solverParams = util.solver.parseParams()
    local linearSolver = util.solver.create(solverParams, gmg)
    if bStokes then
        linearSolver:set_convergence_check(ConvCheck(10000, 1.0e-20, 1.0e-10, true))
    else
        linearSolver:set_convergence_check(ConvCheck(10000, 1.0e-20, 1.0e-8, true))
    end

    local newtonSolver = NewtonSolver()
    newtonSolver:set_linear_solver(linearSolver)
    newtonSolver:set_convergence_check(ConvCheck(500, 1.0e-20, 1.0e-8, true))
    newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
    return newtonSolver
end

function ComputeNonLinearSolution(u, domainDisc, solver)
    util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
    AdjustMeanValue(u, "p")
end

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

print("UG4 nephron lumen Navier-Stokes setup")
print("  grid                  = " .. gridName)
print("  discretization        = " .. discType)
print("  refinements           = " .. numRefs)
print("  MPI processes         = " .. NumProcs())
print("  SI kinematic viscosity= " .. kinematic_viscosity_si .. " m^2/s")
print("  mesh viscosity coeff. = " .. kinematic_viscosity)
print("  reference flow [m3/s] = " .. inflow_ultrafiltrate)
print("  inlet speed [m/s]     = " .. inlet_speed)

local dom = CreateDomain()
local approxSpace = CreateApproxSpace(dom)
local domainDisc = CreateDomainDisc(approxSpace)
local solver = CreateSolver(approxSpace)

print(solver:config_string())

local solution = GridFunction(approxSpace)
solution:set(0.0)
ComputeNonLinearSolution(solution, domainDisc, solver)

local fctCmp = approxSpace:names()
local velCmp = {}
for d = 1, #fctCmp - 1 do
    velCmp[d] = fctCmp[d]
end

local vtkWriter = VTKOutput()
vtkWriter:select(velCmp, "velocity")
vtkWriter:select("p", "pressure")
vtkWriter:print("Results/model_0_lumen_flow", solution)

print("Steady lumen flow solved.")
