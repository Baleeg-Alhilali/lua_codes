-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

--[[
Step 1: lumen velocity with apical water leakage
================================================

This script solves only the fluid problem requested for the first modelling
step:

  * steady incompressible Navier--Stokes in Lumen, using its creeping-flow
    limit for the proximal-tubule Reynolds number;
  * a parabolic ultrafiltrate velocity at Inlet;
  * prescribed water leakage out through Apical;
  * the natural zero-traction / no-stress condition at Outlet.

There is deliberately no Darcy equation and no velocity field in Membrane or
Inter in this file.  The previous coupled version is preserved as
velocity_field_coupled.lua.

The mesh coordinates and velocities use micrometres and seconds.

Literature defaults
-------------------
Human S1 proximal-tubule mean velocity:
  48.15 mm/min = 802.5 um/s, rounded here to 800 um/s.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10117878/

Rabbit proximal-tubule water reabsorption:
  1.9 nL/(mm min).  With a 41.5 um lumen diameter this is approximately
  0.24 um/s through the wall, rounded here to 0.25 um/s.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC436574/

Default parallel run:
  mpirun -np 8 ugshell -ex code/flow/velocity_field.lua

Useful overrides:
  -numRefs 2
  -numPreRefs 0
  -meanInletVelocity 800
  -leakageVelocity 0.25
  -lumenRadius 1.0
  -pressureStabilization 1e-6
  -solverReduction 1e-12
  -inletProjectionCorrection 1.333333  (optional override at numRefs=2)
  -grid ProMeshFiles/runs/3rd_trace_third_boxref5.ugx
  -swc ProMeshFiles/runs/3rd_trace_third_boxref5_outer_radius_1p2.swc
]]

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")

-- The fine-grid iteration is distributed. SuperLU is used only for the small
-- level-0 coarse problem inside geometric multigrid.
AssertPluginsLoaded({"NavierStokes", "SuperLU6", "neuro_collection"})

local totalStart = GetClockS()

--------------------------------------------------------------------------------
-- Configuration
--------------------------------------------------------------------------------

local dim = 3
local gridName = util.GetParam(
    "-grid",
    "ProMeshFiles/runs/3rd_trace_third_boxref5.ugx"
)
local centerlineName = util.GetParam(
    "-swc",
    "ProMeshFiles/runs/3rd_trace_third_boxref5_outer_radius_1p2.swc"
)

-- Two projector refinements are the default for this distributed version.
local numRefs = math.floor(util.GetParamNumber("-numRefs", 2))
local numPreRefs = math.floor(util.GetParamNumber("-numPreRefs", 0))
local meanInletVelocity = util.GetParamNumber("-meanInletVelocity", 800.0) -- um/s
local leakageVelocity = util.GetParamNumber("-leakageVelocity", 0.25)       -- um/s
local lumenRadius = util.GetParamNumber("-lumenRadius", 1.0)                -- um
local pressureStabilization = util.GetParamNumber("-pressureStabilization", 1.0e-6)
local solverReduction = util.GetParamNumber("-solverReduction", 1.0e-12)

-- NavierStokesInflow only constrains the nodes carried by the Inlet subset.
-- On this swept cap, the remaining rim nodes belong to Apical, so the retained
-- P1 flux fraction is (1 - 1/2^level).  Compensate for that trace effect.
local defaultProjectionCorrection = 1.0
if numRefs > 0 then
    local levelFactor = math.pow(2.0, numRefs)
    defaultProjectionCorrection = levelFactor / (levelFactor - 1.0)
end
local inletProjectionCorrection = util.GetParamNumber(
    "-inletProjectionCorrection",
    defaultProjectionCorrection
)

-- Constant-density incompressible form divided by density.  The pressure
-- unknown is therefore kinematic pressure p/rho; rho = 1 keeps the momentum
-- and continuity blocks consistently scaled in micrometre-second units.
local waterDensity = 1.0
local kinematicViscosity = 6.96e5 -- um^2/s

assert(numRefs >= 0, "numRefs must be non-negative")
assert(numRefs >= numPreRefs, "numRefs must be >= numPreRefs")
assert(meanInletVelocity > 0.0, "meanInletVelocity must be positive")
assert(leakageVelocity >= 0.0, "leakageVelocity must be non-negative")
assert(lumenRadius > 0.0, "lumenRadius must be positive")
assert(pressureStabilization > 0.0, "pressureStabilization must be positive")
assert(solverReduction > 0.0, "solverReduction must be positive")
assert(inletProjectionCorrection > 0.0, "inletProjectionCorrection must be positive")

if numRefs < 2 then
    print("WARNING: -numRefs " .. numRefs
        .. " is a validation override; use the default -numRefs 2 for production.")
end

--------------------------------------------------------------------------------
-- Centerline geometry for the curved inlet and wall conditions
--------------------------------------------------------------------------------

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

local inletCenter = centerline[2]
local inletDirection = normalizedDirection(centerline[1], centerline[2])
local activeStartS = centerline[2].s
local activeEndS = centerline[#centerline].s
local leakageRampLength = 2.0 * lumenRadius
local inletPeakVelocity = 2.0 * meanInletVelocity

local function nearestCenterlinePoint(x, y, z)
    local bestX = centerline[1][1]
    local bestY = centerline[1][2]
    local bestZ = centerline[1][3]
    local bestS = 0.0
    local bestR2 = math.huge

    for i = 1, #centerline - 1 do
        local a, b = centerline[i], centerline[i + 1]
        local sx, sy, sz = b[1] - a[1], b[2] - a[2], b[3] - a[3]
        local segmentLength2 = sx * sx + sy * sy + sz * sz
        if segmentLength2 > 0.0 then
            local alpha = ((x - a[1]) * sx + (y - a[2]) * sy + (z - a[3]) * sz)
                / segmentLength2
            alpha = math.max(0.0, math.min(1.0, alpha))
            local qx = a[1] + alpha * sx
            local qy = a[2] + alpha * sy
            local qz = a[3] + alpha * sz
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
    local dx = x - inletCenter[1]
    local dy = y - inletCenter[2]
    local dz = z - inletCenter[3]
    local axial = dx * inletDirection[1] + dy * inletDirection[2] + dz * inletDirection[3]
    local radius2 = math.max(0.0, dx * dx + dy * dy + dz * dz - axial * axial)
    return math.max(0.0, 1.0 - radius2 / (lumenRadius * lumenRadius))
end

function UnitInletAxialVelocity(x, y, z, t, si)
    return inletDirection[1], inletDirection[2], inletDirection[3]
end

function UnitInletParabolicVelocity(x, y, z, t, si)
    local shape = inletShape(x, y, z)
    return shape * inletDirection[1], shape * inletDirection[2], shape * inletDirection[3]
end

-- Scalar forms are used to integrate the P1 representation on the actual cap.
function UnitInletParabolicX(x, y, z, t, si)
    return inletShape(x, y, z) * inletDirection[1]
end

function UnitInletParabolicY(x, y, z, t, si)
    return inletShape(x, y, z) * inletDirection[2]
end

function UnitInletParabolicZ(x, y, z, t, si)
    return inletShape(x, y, z) * inletDirection[3]
end

function UltrafiltrateInletVelocity(x, y, z, t, si)
    local speed = inletPeakVelocity * inletShape(x, y, z)
    return speed * inletDirection[1], speed * inletDirection[2], speed * inletDirection[3]
end

-- Positive radial velocity leaves Lumen through Apical.  Leakage is smoothly
-- suppressed at both cap rims so it does not conflict with the inlet trace.
function ApicalLeakageVelocity(x, y, z, t, si)
    local qx, qy, qz, s, radius2 = nearestCenterlinePoint(x, y, z)
    if radius2 < 1.0e-24 then return 0.0, 0.0, 0.0 end
    local scale = leakageVelocity * leakageRamp(s) / math.sqrt(radius2)
    return scale * (x - qx), scale * (y - qy), scale * (z - qz)
end

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

-- The stabilized equal-order discretization stores (u,v,w,p) together at
-- every vertex, so a four-component CPU block is the valid 3D layout.
InitUG(dim, AlgebraType("CPU", 4))
print("MPI processes: " .. NumProcs())

-- Register projSH before loading so LoadDomain also restores the UGX-embedded
-- ProjectionHandler.  This mesh contains NeuriteProjectors on Basolateral and
-- Apical; they keep both curved tube surfaces on their spline geometry when
-- new vertices are created.
--
-- Register the additional handler so LoadDomain takes the UGX path that reads
-- its serialized ProjectionHandler.  The local UG4 loader broadcasts that
-- handler after installation, giving every MPI rank the same Basolateral and
-- Apical neurite projectors before distributed refinement begins.
local dom = Domain()
dom:create_additional_subset_handler("projSH")
print("Loading and broadcasting the embedded neurite refinement projectors ...")
LoadDomain(dom, gridName)
assert(
    util.CheckSubsets(
        dom,
        {"Lumen", "Membrane", "Basolateral", "Apical", "Inter",
         "OuterWall", "Inlet", "Outlet"}
    ),
    "The domain is missing one or more required subsets"
)
print("The UGX ProjectionHandler was broadcast to every MPI rank; both "
    .. "refinement passes use the same neurite projectors.")

-- Pre-refinement, if requested, must also use the loaded projector.
if numPreRefs > 0 then
    local preRefiner = GlobalDomainRefiner(dom)
    for i = 1, numPreRefs do
        preRefiner:refine()
    end
    delete(preRefiner)
end
balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

local nsSpace = ApproximationSpace(dom)
local lumenSupport = {"Lumen", "Apical", "Inlet", "Outlet"}
nsSpace:add_fct({"u", "v", "w"}, "Lagrange", 1, lumenSupport)
nsSpace:add_fct("p_lumen", "Lagrange", 1, "Lumen,Apical,Inlet,Outlet")
nsSpace:init_levels()
nsSpace:init_top_surface()
nsSpace:print_statistic()

-- Normalize the parabola over the projected inlet geometry and its discrete P1
-- trace, then apply the subset-trace correction calculated above.
local inletIntegrationGridFunction = GridFunction(nsSpace)
local signedInletArea = IntegralNormalComponentOnManifold(
    "UnitInletAxialVelocity", inletIntegrationGridFunction,
    "Inlet", "Lumen", 0.0, 5
)
Interpolate("UnitInletParabolicX", inletIntegrationGridFunction, "u", 0.0)
Interpolate("UnitInletParabolicY", inletIntegrationGridFunction, "v", 0.0)
Interpolate("UnitInletParabolicZ", inletIntegrationGridFunction, "w", 0.0)
local discreteUnitProfile = GridFunctionVectorData(
    inletIntegrationGridFunction, "u,v,w"
)
local signedUnitProfileFlow = IntegralNormalComponentOnManifold(
    discreteUnitProfile, inletIntegrationGridFunction,
    "Inlet", "Lumen", 0.0, 5
)
assert(math.abs(signedUnitProfileFlow) > 1.0e-14, "Degenerate inlet profile")
inletPeakVelocity = meanInletVelocity * signedInletArea / signedUnitProfileFlow
    * inletProjectionCorrection

print(string.format("Inlet area: %.9e um^2", math.abs(signedInletArea)))
print(string.format("Parabolic peak: %.9e um/s", inletPeakVelocity))

--------------------------------------------------------------------------------
-- Navier--Stokes equation and boundary conditions
--------------------------------------------------------------------------------

local nsEq = NavierStokes("u,v,w,p_lumen", "Lumen", "fe")
nsEq:set_density(waterDensity)
nsEq:set_kinematic_viscosity(kinematicViscosity)
nsEq:set_stokes(true)
nsEq:set_laplace(false) -- symmetric-gradient viscous stress
nsEq:set_exact_jacobian(true)
nsEq:set_quad_order(5)
nsEq:set_stabilization(pressureStabilization)

local inletBC = NavierStokesInflow(nsEq)
inletBC:add("UltrafiltrateInletVelocity", "Inlet")

local leakageBC = NavierStokesInflow(nsEq)
leakageBC:add("ApicalLeakageVelocity", "Apical")

local nsDomainDisc = DomainDiscretization(nsSpace)
nsDomainDisc:add(nsEq)
nsDomainDisc:add(inletBC)
nsDomainDisc:add(leakageBC)

-- Outlet is deliberately left unconstrained.  In the FE weak form this is
-- the natural zero-traction / no-stress condition.

--------------------------------------------------------------------------------
-- Distributed iterative solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = nsSpace

-- Overlapping ILU is used on each MPI partition.  Unlike the sequential Vanka
-- correction, this smoother also operates on parallel refinement levels.
-- Pressure stabilization is required for the equal-order P1/P1 pair and also
-- supplies the small pressure-block diagonal needed by ILU.
local parallelSmoother = {
    type = "ilu",
    damping = 0.9,
    overlap = true,
    consistentInterfaces = true,
    -- Pressure stabilization may be below ILU's generic 1e-8 pivot cutoff.
    inversionEps = 1.0e-12
}

local linearSolverDesc = {
    type = "bicgstab",
    precond = {
        type = "gmg",
        approxSpace = nsSpace,
        baseLevel = 0,
        baseSolver = "superlu",
        gatheredBaseSolverIfAmbiguous = true,
        cycle = "F",
        preSmooth = 1,
        postSmooth = 1,
        rap = true,
        smoother = parallelSmoother,
        transfer = "std"
    },
    convCheck = {
        type = "standard",
        iterations = 200,
        absolute = 1.0e-8,
        reduction = solverReduction,
        verbose = true
    }
}

print("Solving low-Reynolds-number Navier--Stokes in Lumen (creeping-flow limit) ...")
local solveStart = GetClockS()
local _, nsSolution = util.solver.SolveLinearProblem(nsDomainDisc, linearSolverDesc)

print(string.format("Lumen solve runtime: %.6f min", (GetClockS() - solveStart) / 60.0))

--------------------------------------------------------------------------------
-- Output and mass-balance diagnostics
--------------------------------------------------------------------------------

local lumenVelocityData = GridFunctionVectorData(nsSolution, "u,v,w")
local lumenPressureData = GridFunctionNumberData(nsSolution, "p_lumen")
-- The explicit name prevents this refined result from being confused with
-- older base-grid files named velocity_field_lumen_t0000.vtu.
local outputBase = "Results/velocity_field_lumen_projected_after_"
    .. numRefs .. "_refinements"

local lumenOut = VTKOutput()
lumenOut:clear_selection()
-- With P1/P1 and CPU blocks of size four, velocity and pressure are both
-- exported at every vertex of the twice-refined top-surface grid.
-- Selecting data objects also works around UG4's empty PPointData metadata
-- when symbolic functions are written using parallel print_subsets().
lumenOut:select(lumenVelocityData, "velocity_um_per_s")
lumenOut:select(lumenPressureData, "kinematic_pressure_um2_per_s2")
print("Writing the top-surface Lumen solution after " .. numRefs
    .. " refinement(s) ...")
lumenOut:print_subsets(outputBase, nsSolution, "Lumen", 0, 0.0)

local inletFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Inlet", "Lumen", 0.0, 5
)
local leakageFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Apical", "Lumen", 0.0, 5
)
local outletFlow = IntegralNormalComponentOnManifold(
    lumenVelocityData, nsSolution, "Outlet", "Lumen", 0.0, 5
)
local balance = inletFlow + leakageFlow + outletFlow

print("Flux convention: positive means outward from Lumen.")
print(string.format("Lumen inlet flux:      %.9e um^3/s", inletFlow))
print(string.format("Lumen apical leakage:  %.9e um^3/s", leakageFlow))
print(string.format("Lumen outlet flux:     %.9e um^3/s", outletFlow))
print(string.format("Lumen balance residual:%.9e um^3/s", balance))
print("Output: " .. outputBase)
print(string.format("Total runtime:         %.6f min", (GetClockS() - totalStart) / 60.0))
