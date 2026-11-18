--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- Solves steady lumen flow followed by one-way-coupled Darcy flow through
-- Membrane and Inter. A prescribed Apical leakage velocity supplies the Darcy
-- influx, while Outlet remains a natural zero-traction boundary.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolve every relative input and output path from the Model directory and
-- load the UG4 libraries and plugins required by both subproblems.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--[[
Steady fluid field in the nephron mesh
======================================

  * incompressible Navier--Stokes in Lumen
  * Darcy flow in Membrane and Inter
  * parabolic ultrafiltrate velocity at Inlet
  * zero-normal-stress (traction-free) condition at Outlet
  * prescribed, literature-based water leakage through Apical into the porous region

The two problems are solved sequentially.  The coupling is one-way: the same
target leakage prescribed out of Lumen is imposed as the inflow flux of the
porous Membrane problem.  Darcy's law and the measured membrane hydraulic
conductivity then predict the effective hydrostatic/osmotic driving pressure.

Units are micrometres, seconds, pascals, and kilograms where applicable.

Literature defaults
-------------------
Human S1 proximal-tubule mean velocity:
  48.15 mm/min = 802.5 um/s, rounded here to 800 um/s.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10117878/

Rabbit proximal-tubule water reabsorption:
  1.9 nL/(mm min).  With a 41.5 um lumen diameter this is approximately
  0.24 um/s through the wall, rounded here to 0.25 um/s.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC436574/

Measured proximal-tubule water hydraulic conductivity:
  2.9--6.3e-5 cm/(s atm); the midpoint is used to estimate membrane K.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC291894/

Runtime inputs (all other values are in the Parameters block below):
  -numRefs 0
  -output Results/velocity_field_coupled
]]

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- This checkout's current plugin exports the UG4 class group "SuperLU" from
-- the dynamically loaded plugin named "SuperLU6".
AssertPluginsLoaded({"NavierStokes", "ConvectionDiffusion", "SuperLU6"})

local totalStart = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------
-- Only -numRefs and -output are command-line inputs. Mesh paths, physical
-- coefficients, discretization choices, and solver choices are set here.

local dim = 3
local gridName = "ProMeshFiles/runs/3rd_trace_third_boxref5.ugx"
local centerlineName =
    "ProMeshFiles/runs/3rd_trace_third_boxref5_outer_radius_1p2.swc"

-- SuperLU is a serial direct solver (gathered when MPI is used).  Keep the
-- direct-solver default at the coarse mesh; increase only if memory permits.
local numRefs = math.floor(util.GetParamNumber("-numRefs", 0))
local numPreRefs = 0
local solveFullNavierStokes = false

local meanInletVelocity = 800.0 -- um/s
local leakageVelocity = 0.25 -- um/s
local lumenRadius = 1.0 -- um
local membraneThickness = 0.2 -- um
-- On the mesh's single-quadrilateral cap, NavierStokesInflowFE transfers
-- two-thirds of the analytically integrated P2 parabola into the solved trace.
local inletProjectionCorrection = 1.5
local pressureStabilization = 1.0e-6
local outputFileName = util.GetParam(
    "-output",
    "Results/velocity_field_coupled",
    "base prefix for lumen, membrane, and inter VTK output")

-- Water at 37 C.
local waterDensity = 993.0                     -- kg/m^3
local dynamicViscosity = 6.91e-4              -- Pa s
local kinematicViscosity = 6.96e5             -- um^2/s

-- Midpoint of 2.9--6.3e-5 cm/(s atm), converted to um/(Pa s).
local membraneHydraulicConductivity = 4.54e-6 -- um/(Pa s)
local membranePermeability =
    membraneHydraulicConductivity * dynamicViscosity * membraneThickness -- um^2

-- This is a modelling parameter because the mesh's Inter region is an
-- effective porous continuum.  It is deliberately much more permeable than
-- the epithelial Membrane.
local interPermeability = 1.0e-3 -- um^2

assert(numRefs >= numPreRefs, "numRefs must be >= numPreRefs")
assert(meanInletVelocity > 0.0, "meanInletVelocity must be positive")
assert(leakageVelocity >= 0.0, "leakageVelocity must be non-negative")
assert(lumenRadius > 0.0, "lumenRadius must be positive")
assert(inletProjectionCorrection > 0.0, "inletProjectionCorrection must be positive")
assert(membranePermeability > 0.0, "membranePermeability must be positive")
assert(interPermeability > 0.0, "interPermeability must be positive")
assert(pressureStabilization > 0.0, "pressureStabilization must be positive")

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- Load and refine the complete multiphysics mesh, then construct the lumen and
-- porous approximation spaces on the shared domain.

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

local dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
     "OuterWall", "Inlet", "Outlet"},
    true
)
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

local nsSpace = ApproximationSpace(dom)
local lumenSupport = {"Lumen", "Apical", "Inlet", "Outlet"}
nsSpace:add_fct({"u", "v", "w"}, "Lagrange", 2, lumenSupport)
nsSpace:add_fct("p_lumen", "Lagrange", 1, "Lumen,Apical,Inlet,Outlet")
nsSpace:init_levels()
nsSpace:init_top_surface()
nsSpace:print_statistic()

local darcySpace = ApproximationSpace(dom)
darcySpace:add_fct(
    "p_darcy",
    "Lagrange",
    1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,Inlet,Outlet"
)
darcySpace:init_levels()
darcySpace:init_top_surface()
darcySpace:print_statistic()

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Define all mesh-dependent boundary callbacks first so the equation blocks
-- below can refer to them by their registered Lua names.

local function loadCenterline(path)
    local points = {}
    local file = assert(io.open(path, "r"), "Cannot open centerline: " .. path)

    for line in file:lines() do
        local payload = line:gsub("#.*$", "")
        local values = {}
        for token in payload:gmatch("%S+") do
            values[#values + 1] = tonumber(token)
        end
        if #values == 7 then
            points[#points + 1] = {values[3], values[4], values[5], s = 0.0}
        end
    end
    file:close()

    assert(#points >= 2, "The centerline must contain at least two SWC points")
    for i = 2, #points do
        local dx = points[i][1] - points[i - 1][1]
        local dy = points[i][2] - points[i - 1][2]
        local dz = points[i][3] - points[i - 1][3]
        points[i].s = points[i - 1].s + math.sqrt(dx * dx + dy * dy + dz * dz)
    end
    return points
end

local centerline = loadCenterline(centerlineName)

local function normalizedDirection(a, b)
    local dx, dy, dz = b[1] - a[1], b[2] - a[2], b[3] - a[3]
    local length = math.sqrt(dx * dx + dy * dy + dz * dz)
    assert(length > 0.0, "Repeated points at the end of the centerline")
    return {dx / length, dy / length, dz / length}
end

-- The importer places the inlet cap at the second SWC sample.  The first
-- segment supplies its inward tangent.
local inletCenter = centerline[2]
local inletDirection = normalizedDirection(centerline[1], centerline[2])
local activeStartS = centerline[2].s
local activeEndS = centerline[#centerline].s
local leakageRampLength = 2.0 * lumenRadius
local inletPeakVelocity = 2.0 * meanInletVelocity
local darcyLeakageScale = 1.0

-- Return the nearest centerline point, its arc coordinate, and squared radius.
local function nearestCenterlinePoint(x, y, z)
    local bestX, bestY, bestZ, bestS = centerline[1][1], centerline[1][2], centerline[1][3], 0.0
    local bestR2 = math.huge

    for i = 1, #centerline - 1 do
        local a, b = centerline[i], centerline[i + 1]
        local sx, sy, sz = b[1] - a[1], b[2] - a[2], b[3] - a[3]
        local segmentLength2 = sx * sx + sy * sy + sz * sz
        if segmentLength2 > 0.0 then
            local alpha = ((x - a[1]) * sx + (y - a[2]) * sy + (z - a[3]) * sz)
                / segmentLength2
            alpha = math.max(0.0, math.min(1.0, alpha))
            local qx, qy, qz = a[1] + alpha * sx, a[2] + alpha * sy, a[3] + alpha * sz
            local rx, ry, rz = x - qx, y - qy, z - qz
            local r2 = rx * rx + ry * ry + rz * rz
            if r2 < bestR2 then
                bestR2 = r2
                bestX, bestY, bestZ = qx, qy, qz
                bestS = a.s + alpha * math.sqrt(segmentLength2)
            end
        end
    end

    return bestX, bestY, bestZ, bestS, bestR2
end

local function leakageRamp(s)
    if leakageVelocity == 0.0 then return 0.0 end
    local distanceFromCap = math.min(s - activeStartS, activeEndS - s)
    return math.max(0.0, math.min(1.0, distanceFromCap / leakageRampLength))
end

local function inletShape(x, y, z)
    local dx, dy, dz = x - inletCenter[1], y - inletCenter[2], z - inletCenter[3]
    local axial = dx * inletDirection[1] + dy * inletDirection[2] + dz * inletDirection[3]
    local radius2 = math.max(0.0, dx * dx + dy * dy + dz * dz - axial * axial)
    return math.max(0.0, 1.0 - radius2 / (lumenRadius * lumenRadius))
end

-- These two unit profiles are integrated below to normalize the parabolic
-- inlet over the actual mesh cap.  This matters here because the coarse cap
-- is a quadrilateral approximation rather than an exact circle.
function UnitInletAxialVelocity(x, y, z, t, si)
    return inletDirection[1], inletDirection[2], inletDirection[3]
end

function UnitInletParabolicVelocity(x, y, z, t, si)
    local shape = inletShape(x, y, z)
    return shape * inletDirection[1], shape * inletDirection[2], shape * inletDirection[3]
end

function UltrafiltrateInletVelocity(x, y, z, t, si)
    local speed = inletPeakVelocity * inletShape(x, y, z)
    return speed * inletDirection[1], speed * inletDirection[2], speed * inletDirection[3]
end

-- Prescribed velocity out of Lumen and into Membrane.  The nearest-centreline
-- vector is the radial surface normal of this projector-generated tube.
function ApicalLeakageVelocity(x, y, z, t, si)
    local qx, qy, qz, s, radius2 = nearestCenterlinePoint(x, y, z)
    if radius2 < 1.0e-24 then return 0.0, 0.0, 0.0 end
    local scale = leakageVelocity * leakageRamp(s) / math.sqrt(radius2)
    return scale * (x - qx), scale * (y - qy), scale * (z - qz)
end

-- UG4's Neumann sign is outward from Membrane.  Water entering Membrane at its
-- inner boundary therefore has a negative sign.  The magnitude uses the
-- identical cap ramp as ApicalLeakageVelocity.
function DarcyLeakageFlux(x, y, z, t, si)
    local _, _, _, s = nearestCenterlinePoint(x, y, z)
    return -darcyLeakageScale * leakageVelocity * leakageRamp(s)
end

function DarcyLeakageVelocity(x, y, z, t, si)
    local vx, vy, vz = ApicalLeakageVelocity(x, y, z, t, si)
    return darcyLeakageScale * vx, darcyLeakageScale * vy, darcyLeakageScale * vz
end

--------------------------------------------------------------------------------
-- Lumen equation
--------------------------------------------------------------------------------

-- Normalize so -integral(u.n)dA = meanInletVelocity * inlet area exactly,
-- independently of the polygonal cap approximation.
local inletIntegrationGridFunction = GridFunction(nsSpace)
local signedInletArea = IntegralNormalComponentOnManifold(
    "UnitInletAxialVelocity", inletIntegrationGridFunction,
    "Inlet", "Lumen", 0.0, 5
)
local signedUnitProfileFlow = IntegralNormalComponentOnManifold(
    "UnitInletParabolicVelocity", inletIntegrationGridFunction,
    "Inlet", "Lumen", 0.0, 5
)
assert(math.abs(signedUnitProfileFlow) > 1.0e-14, "Degenerate inlet profile")
inletPeakVelocity = meanInletVelocity * signedInletArea / signedUnitProfileFlow
    * inletProjectionCorrection
print(string.format(
    "Inlet area %.9e um^2; mesh-corrected parabolic peak %.9e um/s",
    math.abs(signedInletArea), inletPeakVelocity
))

local nsEq = NavierStokes("u,v,w,p_lumen", "Lumen", "fe")
nsEq:set_density(waterDensity)
nsEq:set_kinematic_viscosity(kinematicViscosity)
-- Start from the creeping-flow member of Navier--Stokes.  At the default
-- parameters Re = U D / nu is only about 0.0023, so this is also an excellent
-- physical approximation and a robust initial iterate for the full equations.
nsEq:set_stokes(true)
nsEq:set_laplace(false) -- symmetric-gradient viscous stress
nsEq:set_exact_jacobian(true)
nsEq:set_quad_order(5)
-- This swept mesh has only one hexahedron across the lumen.  A small pressure
-- Laplacian suppresses the resulting unresolved cross-sectional pressure mode.
nsEq:set_stabilization(pressureStabilization)

local inletBC = NavierStokesInflow(nsEq)
inletBC:add("UltrafiltrateInletVelocity", "Inlet")

-- This boundary helper prescribes all three velocity components.  It is used
-- here as a permeable-wall condition; positive radial velocity leaves Lumen.
local leakageBC = NavierStokesInflow(nsEq)
leakageBC:add("ApicalLeakageVelocity", "Apical")

local nsDomainDisc = DomainDiscretization(nsSpace)
nsDomainDisc:add(nsEq)
nsDomainDisc:add(inletBC)
nsDomainDisc:add(leakageBC)
-- For the FE weak form, leaving all outlet fields unconstrained supplies the
-- natural zero-traction / no-stress condition.  UG4's explicit outflow helper
-- is available only for its finite-volume Navier--Stokes discretizations.

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- Solve the lumen and porous problems sequentially. SuperLU is retained in
-- this legacy one-way model; the monolithic production scripts use GMG.

util.solver.defaults.approxSpace = nsSpace
local nsLinearSolverDesc = {
    type = "superlu" -- alternative: an iterative descriptor for large systems
}
local nsSolverDesc = {
    type = "newton", -- alternative for Stokes-only runs: a linear solver
    lineSearch = {
        type = "standard",
        maxSteps = 10,
        lambdaStart = 1.0,
        lambdaReduce = 0.5,
        acceptBest = true,
        checkAll = false
    },
    convCheck = {
        type = "standard",
        iterations = 20,
        absolute = 1e-14,
        reduction = 1.0e-10,
        verbose = true
    },
    linSolver = "superlu" -- alternative: an iterative descriptor for large systems
}

print("Solving the Stokes initial field in Lumen ...")
local nsStart = GetClockS()
local _, nsSolution = util.solver.SolveLinearProblem(nsDomainDisc, nsLinearSolverDesc)

if solveFullNavierStokes then
    -- Add the convective term and converge the full steady Navier--Stokes
    -- system from the Stokes field.  At the default Re this correction is
    -- negligible but can be requested explicitly for comparison.
    nsEq:set_stokes(false)
    local nsSolver = util.solver.CreateSolver(nsSolverDesc)
    nsSolver:init(AssembledOperator(nsDomainDisc))
    nsSolver:prepare(nsSolution)
    print("Applying the full Navier--Stokes correction in Lumen ...")
    assert(nsSolver:apply(nsSolution), "Navier--Stokes solve failed")
else
    print("Using the creeping-flow Navier--Stokes limit configured in Parameters.")
end
print(string.format("Navier--Stokes runtime: %.6f min", (GetClockS() - nsStart) / 60.0))

local lumenVelocityData = GridFunctionVectorData(nsSolution, "u,v,w")
local realizedLumenLeakFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Apical", "Lumen", 0.0, 5
)
local nominalLumenLeakFlow = IntegralNormalComponentOnManifold(
    "ApicalLeakageVelocity", nsSolution, "Apical", "Lumen", 0.0, 5
)
if math.abs(nominalLumenLeakFlow) > 1.0e-14 then
    darcyLeakageScale = realizedLumenLeakFlow / nominalLumenLeakFlow
end
print(string.format("Darcy leakage projection scale: %.9e", darcyLeakageScale))

--------------------------------------------------------------------------------
-- Darcy discretization and boundary conditions
--------------------------------------------------------------------------------

-- Two constant-coefficient element discs share p_darcy at Basolateral.  Their
-- assembled contributions enforce the usual pressure and normal-flux matching
-- across the material interface without a fragile Lua coefficient callback.
local membraneDarcyEq = ConvectionDiffusionFV1("p_darcy", "Membrane")
membraneDarcyEq:set_mass_scale(0.0)
membraneDarcyEq:set_diffusion(membranePermeability / dynamicViscosity)

local interDarcyEq = ConvectionDiffusionFV1("p_darcy", "Inter")
interDarcyEq:set_mass_scale(0.0)
interDarcyEq:set_diffusion(interPermeability / dynamicViscosity)

local membraneDarcyVelocity = DarcyVelocityLinker()
membraneDarcyVelocity:set_permeability(membranePermeability)
membraneDarcyVelocity:set_viscosity(dynamicViscosity)
membraneDarcyVelocity:set_density(waterDensity)
membraneDarcyVelocity:set_gravity(ConstUserVector(0.0))
membraneDarcyVelocity:set_pressure_gradient(membraneDarcyEq:gradient())

local interDarcyVelocity = DarcyVelocityLinker()
interDarcyVelocity:set_permeability(interPermeability)
interDarcyVelocity:set_viscosity(dynamicViscosity)
interDarcyVelocity:set_density(waterDensity)
interDarcyVelocity:set_gravity(ConstUserVector(0.0))
interDarcyVelocity:set_pressure_gradient(interDarcyEq:gradient())

-- Prescribe the water entering Membrane through Apical.  The resulting Apical
-- pressure is predicted from the membrane resistance, rather than prescribed.
local darcyLeakageBC = NeumannBoundary("p_darcy", "fv1")
darcyLeakageBC:add("DarcyLeakageFlux", "Apical", "Membrane")

-- The outer interstitial box is the zero effective-pressure reference and the
-- exit for reabsorbed water.  All other porous boundaries are no-flux.
local darcyReference = DirichletBoundary()
darcyReference:add(0.0, "p_darcy", "OuterWall")

local darcyDomainDisc = DomainDiscretization(darcySpace)
darcyDomainDisc:add(membraneDarcyEq)
darcyDomainDisc:add(interDarcyEq)
darcyDomainDisc:add(darcyLeakageBC)
darcyDomainDisc:add(darcyReference)

util.solver.defaults.approxSpace = darcySpace
local darcySolverDesc = {
    type = "superlu" -- alternative: an iterative descriptor for large systems
}
print("Solving Darcy flow in Membrane and Inter ...")
local darcyStart = GetClockS()
local _, darcySolution = util.solver.SolveLinearProblem(darcyDomainDisc, darcySolverDesc)
print(string.format("Darcy runtime: %.6f min", (GetClockS() - darcyStart) / 60.0))

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- Both subproblems are steady. The linear solves create their zero initial
-- states internally; the optional nonlinear correction starts from Stokes.

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- Write separate lumen, membrane, and inter files from one requested prefix,
-- then report the fluxes used to check conservation and one-way coupling.

local lumenOut = VTKOutput()
lumenOut:clear_selection()
-- Write the P2 velocity at cell centres.  A nodal-only export would omit the
-- P2 edge/interior values and hide most of the axial flow on this one-cell-
-- across lumen mesh.
lumenOut:select_element(lumenVelocityData, "velocity_um_per_s")
lumenOut:select_nodal("p_lumen", "pressure_lumen")
lumenOut:print_subsets(outputFileName .. "_lumen", nsSolution, "Lumen", 0, 0.0)

local membraneOut = VTKOutput()
membraneOut:clear_selection()
membraneOut:select_element(membraneDarcyVelocity, "darcy_velocity_um_per_s")
membraneOut:select_nodal("p_darcy", "effective_pressure_Pa")
membraneOut:print_subsets(
    outputFileName .. "_membrane",
    darcySolution,
    "Membrane",
    0,
    0.0
)

local interOut = VTKOutput()
interOut:clear_selection()
interOut:select_element(interDarcyVelocity, "darcy_velocity_um_per_s")
interOut:select_nodal("p_darcy", "effective_pressure_Pa")
interOut:print_subsets(
    outputFileName .. "_inter",
    darcySolution,
    "Inter",
    0,
    0.0
)

local inletFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Inlet", "Lumen", 0.0, 5
)
local lumenLeakFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Apical", "Lumen", 0.0, 5
)
local lumenOutletFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Outlet", "Lumen", 0.0, 5
)
local prescribedDarcyLeakFlow = IntegralNormalComponentOnManifold(
    "DarcyLeakageVelocity", darcySolution, "Apical", "Membrane", 0.0, 5
)

print("Flux convention: positive means outward from the named inner subset.")
print(string.format("Lumen inlet flux:       %.9e um^3/s", inletFlow))
print(string.format("Lumen apical leakage:   %.9e um^3/s", lumenLeakFlow))
print(string.format("Lumen outlet flux:      %.9e um^3/s", lumenOutletFlow))
print(string.format("Membrane imposed influx:%.9e um^3/s", prescribedDarcyLeakFlow))
print(string.format("Total runtime:          %.6f min", (GetClockS() - totalStart) / 60.0))
