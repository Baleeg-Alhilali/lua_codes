--------------------------------------------------------------------------------
-- Model summary
--------------------------------------------------------------------------------
-- This script solves steady Stokes--Darcy flow and then a one-way-coupled
-- sodium transport problem driven by the converged velocity field:
--
--   * Stokes flow is assembled only in the Lumen volume.
--   * Darcy flow is assembled in the Membrane and Inter volumes.
--   * The Apical surface couples both models through pressure continuity,
--     normal mass conservation, and the Beavers--Joseph--Saffman tangential
--     condition. A geometry-local tangential resistance represents the
--     stationary epithelial wall without selecting a global Cartesian
--     velocity component; the normal component remains free for leakage.
--   * One sodium concentration is transported by the Stokes velocity in the
--     lumen and by the Darcy velocities in membrane and interstitium.
--   * Passive concentration leaks connect Lumen/Membrane at Apical and
--     Membrane/Inter at Basolateral. No electrical migration is included.
--   * p_stokes is the physical kinematic pressure. p_darcy_scaled stores
--     p_D / darcyPressureScale for matrix conditioning; output multiplies it
--     back by darcyPressureScale. The scaling does not change the physics.
--
-- PressureTraceConstraint eliminates the duplicate Darcy trace on Apical and
-- imposes p_stokes = darcyPressureScale * p_darcy_scaled. It also merges the
-- Stokes and Darcy conservation rows, which enforces normal-flux continuity
-- without introducing a zero-diagonal Lagrange-multiplier field.
--
-- This is the preserved formulation that solved numRefs=2 and produced the
-- pressure_trace_darcy_fe_ref2 output. The known refined-grid Darcy-pressure
-- oscillation is intentionally not hidden or filtered in this script.
--
-- This independent variant uses K_membrane = K_inter = 1e-12 m^2 and a
-- natural no-flow condition on OuterWall. The outer pressure parameter below
-- is retained only for the optional fixed-pressure variant.

--------------------------------------------------------------------------------
-- Dependencies
--------------------------------------------------------------------------------
-- Resolve the Model directory from this file instead of relying on the shell's
-- current directory. Load only the UG4 utility scripts and plugins used by the
-- mesh hierarchy, finite-volume Stokes operator, finite-element Darcy operator,
-- NeuriteProjector, and GMG flow solve with a coarse-grid SuperLU correction.
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

AssertPluginsLoaded({
    "NavierStokes", "ConvectionDiffusion", "neuro_collection", "SuperLU6"
})

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------
-- The mesh, physical model, units, boundary data, algebraic scaling, and
-- solver tolerances live here so a run can be reproduced from this file alone.

dim = 3

numRefs = util.GetParamNumber("-numRefs", 1, "number of uniform refinements")
flowOnly = util.HasParamOption("-flowOnly")
vmax = 3.0                          -- maximum inlet velocity [mm/s]
naInlet = util.GetParamNumber("-naInlet", 1.0, "inlet Na concentration [mM]")
naOuter = util.GetParamNumber("-naOuter", 0.0, "OuterWall Na concentration [mM]")
naDiffusion = util.GetParamNumber(
    "-naDiffusion", 1.33e-3, "Na diffusion coefficient [mm^2/s]"
)
naBasolateralPermeability = util.GetParamNumber(
    "-naBasolateralPermeability", 5.0e-4,
    "Basolateral Na permeability [mm/s]"
)

-- Geometry and physical parameters for runs/third_trace_new_build.ugx
-- (October 5 mesh): the Inlet contains one center vertex and a 10-vertex
-- circular wall ring.
gridS = "third_trace_new_build"
gridName = util.GetParam(
    "-grid", "runs/" .. gridS .. ".ugx",
    "input nephron UGX file"
)
-- Exact center vertex and designed wall-ring radius of the October 5 mesh.
-- Its cap lies in z=0.2108 and the first centerline segment enters along +z.
-- With this radius the parabolic profile is exactly zero on the wall ring.
inletRadius = 5.0e-4                -- [mm]
inletCenterX = 0.3248241776875      -- [mm]
inletCenterY = 0.46847242884374996  -- [mm]
inletCenterZ = 0.2108               -- [mm]
inletNormalX = 0.0
inletNormalY = 0.0
inletNormalZ = 1.0
viscosity = 1.0                   -- kinematic viscosity [mm^2/s]
membranePermeabilityM2 = 1.0e-14  -- membrane Darcy permeability [m^2]
interPermeabilityM2 = 1.0e-12     -- interstitial Darcy permeability [m^2]
bjsInterfacePermeabilityM2 = 1.0e-12 -- effective BJS slip permeability [m^2]
bjsAlpha = 0.2                    -- dimensionless BJS coefficient
density = 1000.0                  -- [kg/m^3]
outerPressurePa = 2.0             -- [Pa], used only if outerWallNoFlux=false
outerWallNoFlux = true            -- natural Darcy condition q_D.n = 0 on OuterWall
darcyPressureScale = 1.0e4        -- algebraic only; output is rescaled physically
normalVelocityStabilization = 1.0e3
                                    -- Apical normal impedance [mm/s]
apicalTangentialResistance = util.GetParamNumber(
    "-apicalTangentialResistance", 1.0e6,
    "geometry-local Apical tangential resistance [mm/s]"
)
gmgBaseLevel = 0                  -- retain the loaded coarse mesh as the GMG base
flowAbsoluteTolerance = 1.0e-7   -- attainable absolute floor; relative reduction is checked separately
flowSolverVariant = util.GetParam(
    "-flowSolverVariant", "fine_superlu",
    "base_flow, curved_flow, bicgstab_no_rap, gmres_rap, gmres_no_rap, or fine_superlu"
)
local outerWallTag = outerWallNoFlux and "_outerwall_noflux" or ""

-- The Na writer creates one .pvtu master containing all three concentration
-- regions; MPI runs use one referenced .vtu piece per rank and region.
vtk_file_name = util.GetParam(
    "-output",
        "Results_final/coupled_transport_flow_only_Na_leak"
        .. "_curved_nephron_projector"
        .. "_stokes_darcy_" .. flowSolverVariant .. "_na_gmres"
        .. outerWallTag
        .. "_ref" .. numRefs,
    "VTK output prefix"
)


InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Mesh and approximation space
--------------------------------------------------------------------------------
-- The UGX contains the coarse mesh, subset layout, NeuriteProjector, and the
-- npSurfParams attachment required to place new surface vertices.

local function file_contains(path, needle)
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

assert(file_contains(gridName, '<projection_handler'),
    "The nephron UGX has no serialized ProjectionHandler: " .. gridName)
assert(file_contains(gridName, '<projector type="neurite"'),
    "The nephron UGX has no serialized NeuriteProjector: " .. gridName)
assert(file_contains(gridName, '<vertex_attachment name="npSurfParams"'),
    "The nephron UGX has no npSurfParams attachment: " .. gridName)

-- The first load registers the custom npSurfParams attachment type. The real
-- load can then restore every per-vertex neurite coordinate together with the
-- serialized curved-nephron projector.
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
    "OuterWall", "MembraneWall", "Inlet", "Outlet", "TripleJunction"
}), "The curved nephron is missing a required flow/transport subset")

-- Refine the loaded hierarchy directly so each new surface vertex is placed
-- by NeuriteProjector. Distribution happens at the requested finest level;
-- this preserves one identical curved hierarchy for geometric multigrid.
balancer.parametersParsed = true
balancer.firstDistLvl = numRefs
balancer.Rebalance(dom)
if numRefs > 0 then
    local projectorRefiner = GlobalDomainRefiner(dom)
    for level = 1, numRefs do
        projectorRefiner:refine()
        balancer.Rebalance(dom)
    end
    delete(projectorRefiner)
end

--------------------------------------------------------------------------------
-- Approximation-space fields
--------------------------------------------------------------------------------
-- All five unknowns belong to one ApproximationSpace and therefore one global
-- matrix. Their subset support is deliberately wider than their physical PDE
-- domains near Apical: UG4 user-data evaluation needs both pressure traces on
-- elements adjacent to the interface. The non-physical extensions are pinned
-- to zero below and do not acquire Stokes or Darcy volume equations.

approxSpace = ApproximationSpace(dom)
local allFlowSubsets =
    "Lumen,Membrane,Inter,Apical,Basolateral,OuterWall,MembraneWall,Inlet,Outlet,TripleJunction"
local stokesVelocitySubsetNames = {
    "Lumen", "Apical", "Inlet", "Outlet", "TripleJunction"
}
local stokesPressureAndInterfaceSubsetNames = {
    "Lumen", "Membrane", "Apical", "Basolateral", "MembraneWall",
    "Inlet", "Outlet", "TripleJunction"
}
approxSpace:add_fct(
    {"u", "v", "w"}, "Lagrange", 1, stokesVelocitySubsetNames
)
-- p_stokes has a pinned membrane-side extension because the inner-boundary
-- pressure constraint is evaluated from both volumes adjacent to Apical.
approxSpace:add_fct(
    {"p_stokes"}, "Lagrange", 1, stokesPressureAndInterfaceSubsetNames
)
approxSpace:add_fct(
    "p_darcy_scaled", "Lagrange", 1, allFlowSubsets
)
approxSpace:init_levels()
approxSpace:init_top_surface()

-- Keep UG4's native component ordering.  The scalar fields live on different
-- subset unions, so a global lexicographic permutation is not defined for all
-- eight components and is not needed by the ILU--Uzawa multigrid smoother.

-- Print UG4's global distributed-domain statistics followed by the algebraic
-- DoF statistics. The first table reports the mesh at every refinement level;
-- the second reports the unknown count contributed by each solution field.
print("=== MESH STATISTICS ===")
print(string.format("Mesh file: %s", gridName))
print(string.format("MPI processes: %d", NumProcs()))
print(dom:domain_info():to_string())
print("=== APPROXIMATION-SPACE / DOF STATISTICS ===")
approxSpace:print_statistic("component")

--------------------------------------------------------------------------------
-- Discretization and boundary conditions
--------------------------------------------------------------------------------
-- Convert SI permeability to mm^2, construct the parabolic inlet profile, and
-- assemble the Stokes, Darcy, and interface equations. Pressures inside UG4 are
-- kinematic pressures [mm^2/s^2]; multiplication by density*1e-6 converts them
-- to Pa for reporting.

local membranePermeability = membranePermeabilityM2 * 1.0e6 -- [mm^2]
local interPermeability = interPermeabilityM2 * 1.0e6       -- [mm^2]
local bjsInterfacePermeability = bjsInterfacePermeabilityM2 * 1.0e6
local membraneDarcyMobility = membranePermeability / viscosity
local interDarcyMobility = interPermeability / viscosity
local scaledMembraneDarcyMobility =
    membraneDarcyMobility * darcyPressureScale
local scaledInterDarcyMobility = interDarcyMobility * darcyPressureScale
local bjsFriction =
    bjsAlpha * viscosity / math.sqrt(bjsInterfacePermeability)
local outerPressurePhysical = outerPressurePa / (density * 1.0e-6)
local outerPressureScaled = outerPressurePhysical / darcyPressureScale

function inletVelocity(x, y, z, t)
    -- Measure radius in the plane perpendicular to the local inlet direction.
    -- This remains valid when a future nephron inlet is not aligned with z.
    local dx = x - inletCenterX
    local dy = y - inletCenterY
    local dz = z - inletCenterZ
    local axial = dx * inletNormalX + dy * inletNormalY + dz * inletNormalZ
    local rx = dx - axial * inletNormalX
    local ry = dy - axial * inletNormalY
    local rz = dz - axial * inletNormalZ
    local radiusSquared = rx * rx + ry * ry + rz * rz
    local profile = math.max(
        0.0, 1.0 - radiusSquared / (inletRadius * inletRadius)
    )
    if radiusSquared >= inletRadius * inletRadius * (1.0 - 1.0e-6) then
        profile = 0.0
    end
    return vmax * profile * inletNormalX,
           vmax * profile * inletNormalY,
           vmax * profile * inletNormalZ
end

-- Stokes equation in Lumen
--------------------------------------------------------------------------------
-- FV1 uses equal-order P1 velocity and pressure. The fields/cor stabilization
-- suppresses the equal-order pressure null mode. set_laplace(false) selects the
-- symmetric-gradient stress tensor needed by the BJS traction condition.

stokesEq = NavierStokesFV1({"u", "v", "w", "p_stokes"}, {"Lumen"})
stokesEq:set_stokes(true)
stokesEq:set_laplace(false) -- symmetric-gradient stress required for BJS
stokesEq:set_kinematic_viscosity(viscosity)
stokesEq:set_exact_jacobian(true)
stokesEq:set_stabilization("fields", "cor")

-- Prescribe the full velocity at Inlet. Outlet receives the standard zero
-- normal-stress condition. The same open-boundary discretization is required
-- on Apical to assemble the boundary-normal flux in Stokes continuity.
stokesInlet = NavierStokesInflow(stokesEq)
stokesInlet:add("inletVelocity", "Inlet")
stokesOpenBoundary = NavierStokesNoNormalStressOutflow(stokesEq)
-- This boundary discretization is also required on Apical: besides its
-- traction terms, it closes the Stokes continuity control volumes with the
-- actual boundary-normal velocity.  WSBC alone supplies only the tangential
-- BJS term and does not assemble that mass flux.
stokesOpenBoundary:add("Outlet,Apical")
stokesRim = NavierStokesWall(stokesEq)
stokesRim:add("TripleJunction")

-- BJS--Saffman tangential law:
--     (sigma_S n)_t = -alpha*nu/sqrt(K_interface) * u_S,t.
-- The normal factor is a Robin impedance used by the validated formulation to
-- damp an otherwise weak interface-normal mode; it does not replace the normal
-- conservation equation imposed by PressureTraceConstraint.
stokesBJS = NavierStokesWSBCFV1(stokesEq)
stokesBJS:add("Apical")
stokesBJS:set_sliding_factor(math.max(
    bjsFriction, apicalTangentialResistance
))
stokesBJS:set_sliding_limit(0.0)
stokesBJS:set_normal_factor(normalVelocityStabilization)

--------------------------------------------------------------------------------
-- Darcy equation through Membrane and Inter with a scaled pressure unknown
--------------------------------------------------------------------------------
-- With q_D = -(K/nu) grad(p_D), substituting
-- p_D = darcyPressureScale*p_darcy_scaled multiplies the diffusion coefficient
-- by the same scale. Separate coefficients represent the lower-permeability
-- membrane and higher-permeability interstitium. Basolateral is an internal
-- Membrane/Inter interface. OuterWall is either impermeable (the present
-- test) or assigned the fixed pressure configured above.

darcyMembraneEq = ConvectionDiffusionFE("p_darcy_scaled", "Membrane")
darcyMembraneEq:set_mass_scale(0.0)
darcyMembraneEq:set_diffusion(scaledMembraneDarcyMobility)

darcyInterEq = ConvectionDiffusionFE("p_darcy_scaled", "Inter")
darcyInterEq:set_mass_scale(0.0)
darcyInterEq:set_diffusion(scaledInterDarcyMobility)

-- With no explicit boundary discretization, the Darcy diffusion operator uses
-- its natural homogeneous Neumann condition q_D.n=0 on OuterWall.  The fixed
-- pressure alternative remains available for a controlled comparison.
if not outerWallNoFlux then
    darcyOuterPressure = DirichletBoundary()
    darcyOuterPressure:add(
        outerPressureScaled, "p_darcy_scaled", "OuterWall"
    )
end

coupledDisc = DomainDiscretization(approxSpace)
coupledDisc:add(stokesEq)
coupledDisc:add(stokesInlet)
coupledDisc:add(stokesOpenBoundary)
coupledDisc:add(stokesRim)
coupledDisc:add(stokesBJS)

-- Normal traction balance supplies Darcy pressure to the Stokes momentum rows.
-- The neglected viscous normal-stress correction leaves sigma_S*n_S=-p_D*n_S.
-- NeumannBoundary uses the outward conservative momentum flux -sigma_S*n_S,
-- so the applied vector traction is +p_D times each Cartesian basis vector.
local darcyPhysicalPressureValue =
    darcyPressureScale * darcyMembraneEq:value()
for _, row in ipairs({
    {"u", util.CreateConstUserVector3d(1.0, 0.0, 0.0)},
    {"v", util.CreateConstUserVector3d(0.0, 1.0, 0.0)},
    {"w", util.CreateConstUserVector3d(0.0, 0.0, 1.0)}
}) do
    local interfaceTraction = NeumannBoundary(row[1], "fv1")
    interfaceTraction:add(
        darcyPhysicalPressureValue * row[2],
        "Apical", "Lumen"
    )
    coupledDisc:add(interfaceTraction)
end

-- Eliminate the duplicate Darcy trace on Apical and on its endpoint rim:
--     p_darcy_scaled = p_stokes / darcyPressureScale.
-- The eliminated Darcy conservation row is added to the Stokes continuity row.
-- Consequently this one operation imposes physical pressure continuity and
-- balances the Stokes and Darcy normal fluxes without double-counting either.
pressureContinuity = PressureTraceConstraint(
    "p_stokes", "p_darcy_scaled", "Apical,TripleJunction"
)
pressureContinuity:set_pressure_scale(darcyPressureScale)

-- This velocity handle is used only for flux diagnostics and VTK output. No
-- additional Darcy Neumann condition is added on Apical because the merged
-- conservation row already supplies the normal-flux coupling.
local stokesVelocity = stokesEq:velocity()

-- Coupled user data must see both fields on the elements adjacent to Apical.
-- Pin only their non-physical extensions; never pin the Apical trace or the
-- TripleJunction rim, where both physical pressures must remain continuous.
inactiveDarcy = DirichletBoundary()
inactiveDarcy:add(
    0.0, "p_darcy_scaled", "Lumen,Inlet,Outlet"
)
inactiveStokes = DirichletBoundary()
inactiveStokes:add(
    0.0, "p_stokes", "Membrane,Basolateral,MembraneWall"
)

coupledDisc:add(darcyMembraneEq)
coupledDisc:add(darcyInterEq)
coupledDisc:add(pressureContinuity)
if not outerWallNoFlux then coupledDisc:add(darcyOuterPressure) end
coupledDisc:add(inactiveDarcy)
coupledDisc:add(inactiveStokes)

-- Projector refinement creates a valid nonconforming hierarchy.  FV1 needs
-- one-sided P1 constraints to replace every hanging-node slave row by its
-- interpolation relation on all levels used by geometric multigrid.
hangingNodeConstraints = OneSideP1Constraints()
coupledDisc:add(hangingNodeConstraints)
--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------
-- The coupled matrix is indefinite because p_stokes is the Stokes
-- incompressibility multiplier.  The base_flow variant reproduces the working
-- flow code's BiCGStab/GMG/Uzawa configuration exactly. The temporary default
-- is the only variant verified on the complete refined system: fine SuperLU.
--
-- SuperLU is restricted to the fixed level-0 system; the refined system
-- remains iterative GMG.

local solve_start_time = GetClockS()
local function outerLinearSolverDescription(space, disc)
    local gmgMainInverse = ILU()
    gmgMainInverse:set_damp(
        flowSolverVariant == "curved_flow" and 0.8 or 1.0
    )
    gmgMainInverse:set_inversion_eps(1.0e-30)
    local gmgPressureInverse = ILU()
    gmgPressureInverse:set_damp(
        flowSolverVariant == "curved_flow" and 0.8 or 1.0
    )
    gmgPressureInverse:set_inversion_eps(1.0e-30)
    local gmgSmoother = UzawaBase({"p_stokes"})
    gmgSmoother:set_forward_iter(gmgMainInverse)
    gmgSmoother:set_schur_iter(gmgPressureInverse)
    gmgSmoother:set_backward_iter(gmgMainInverse)

    local gmgSchurMassEq = ConvectionDiffusionFV1("p_stokes", "Lumen")
    gmgSchurMassEq:set_mass_scale(1.0)
    gmgSchurMassEq:set_diffusion(0.0)
    local gmgSchurMassDisc = DomainDiscretization(space)
    gmgSchurMassDisc:add(gmgSchurMassEq)
    local gmgSchurMassOp = AssembledLinearOperator(gmgSchurMassDisc)
    gmgSmoother:set_schur_operator_update(gmgSchurMassOp, 1.0 / viscosity)

    local useGMRES = flowSolverVariant == "gmres_rap"
        or flowSolverVariant == "gmres_no_rap"
    local useRAP = flowSolverVariant == "base_flow"
        or flowSolverVariant == "curved_flow"
        or flowSolverVariant == "gmres_rap"
    assert(flowSolverVariant == "base_flow"
        or flowSolverVariant == "curved_flow"
        or flowSolverVariant == "bicgstab_no_rap"
        or flowSolverVariant == "gmres_rap"
        or flowSolverVariant == "gmres_no_rap",
        "Unknown flowSolverVariant: " .. flowSolverVariant)

    return {
        type = useGMRES and "gmres" or "bicgstab",
        precond = {
            type = "gmg",
            approxSpace = space,
            discretization = disc,
            baseLevel = gmgBaseLevel,
            baseSolver = "superlu",
            smoother = gmgSmoother,
            cycle = "V",
            preSmooth = 2,
            postSmooth = 2,
            rap = useRAP,
            gatheredBaseSolverIfAmbiguous = true,
            transfer = {
                type = "std",
                restrictionDamp = 1.0,
                prolongationDamp = 1.0,
                enableP1LagrangeOptimization = true
            }
        },
        convCheck = {
            type = "standard",
            -- These are the exact working-flow linear convergence settings.
            iterations = 400,
            absolute = 1.0e-8,
            reduction = flowSolverVariant == "curved_flow" and 1.0e-4
                or 1.0e-8,
            verbose = true
        }
    }
end

local function solveLinearSubproblem(space, disc, label)
    -- Assemble once because the present Stokes--Darcy problem is linear and
    -- steady. The solution starts at zero after essential BC adjustment.
    util.solver.defaults.approxSpace = space
    local matrix = AssembledLinearOperator(disc)
    local rhs = GridFunction(space)
    local solution = GridFunction(space)
    solution:set(0.0)
    disc:adjust_solution(solution)
    disc:assemble_linear(matrix, rhs)
    local componentNames = {
        "u", "v", "w", "p_stokes", "p_darcy_scaled"
    }
    -- Recompute r=b-Ax from the assembled operator. This guards against a
    -- recursive Krylov residual appearing converged while one physical block
    -- remains inaccurate.
    local function trueDefectAndComponentNorms()
        solution:set_consistent_storage_type()
        local defect = rhs:clone()
        VecAssign(defect, rhs)
        matrix:apply_sub(defect, solution)
        local norms = {}
        for selected = 1, #componentNames do
            local factors = {}
            for component = 1, #componentNames do factors[component] = 0.0 end
            factors[selected] = 1.0
            local componentDefect = defect:clone()
            ScaleGF(componentDefect, defect, factors)
            norms[selected] = VecNorm(componentDefect)
        end
        return defect, norms
    end
    local initialDefectField, initialDefects =
        trueDefectAndComponentNorms()
    local initialGlobalScaleSquared = 0.0
    for component = 1, #componentNames do
        initialGlobalScaleSquared = initialGlobalScaleSquared
            + initialDefects[component] * initialDefects[component]
    end
    local initialGlobalScale = math.sqrt(initialGlobalScaleSquared)
    local solverDesc
    if flowSolverVariant == "fine_superlu" then
        -- Diagnostic only: factor the complete refined matrix, not level 0.
        solverDesc = {type = "superlu"}
    else
        solverDesc = outerLinearSolverDescription(space, disc)
    end
    local solver = util.solver.CreateLinearSolver(solverDesc)
    print("=== " .. label .. " " .. flowSolverVariant
        .. " convergence ===")
    print("=== CONVERGENCE CRITERIA ===")
    print("Configured linear iterations:          400")
    print(string.format(
        "Verified component defect:             <= %.6e",
        flowAbsoluteTolerance
    ))
    print("GMG variants use the working flow solver tolerances.")
    print("Maximum restarted solver cycles:       8")
    print(string.format(
        "Each recomputed component defect:      <= %.6e",
        flowAbsoluteTolerance
    ))
    print("A cycle passes only when every component is PASS.")
    print("=== INITIAL TRUE COMPONENT DEFECTS ===")
    for component = 1, #componentNames do
        print(string.format(
            "INITIAL DEFECT %-14s abs=%.6e",
            componentNames[component], initialDefects[component]
        ))
    end
    print(string.format(
        "INITIAL COMBINED DEFECT              abs=%.6e",
        initialGlobalScale
    ))
    solver:init(matrix, solution)
    local converged = false
    local finalDefectField = initialDefectField
    for restart = 1, 8 do
        print(string.format(
            "=== Solver cycle %d/8 ===", restart
        ))
        solver:apply(solution, rhs)
        local currentDefects
        finalDefectField, currentDefects =
            trueDefectAndComponentNorms()
        local currentGlobalSquared = 0.0
        for component = 1, #componentNames do
            currentGlobalSquared = currentGlobalSquared
                + currentDefects[component] * currentDefects[component]
        end
        local currentGlobal = math.sqrt(currentGlobalSquared)
        local globalRelative = currentGlobal / initialGlobalScale
        print(string.format(
            "TRUE COMBINED DEFECT initial=%.6e current=%.6e "
                .. "relative=%.6e absolute_target=%.6e "
                .. "relative_target=1.000000e-08",
            initialGlobalScale, currentGlobal, globalRelative,
            flowAbsoluteTolerance
        ))
        converged = true
        for component = 1, #componentNames do
            -- Darcy has no independent right-hand side: its initial component
            -- defect is exactly zero and is generated only after Stokes--Darcy
            -- coupling becomes active. Use the global initial defect as the
            -- reference floor so this component receives the same meaningful
            -- relative tolerance as every non-zero initial component.
            local reference = math.max(
                initialDefects[component], initialGlobalScale
            )
            local relative = currentDefects[component] / reference
            local passed =
                currentDefects[component] <= flowAbsoluteTolerance
            print(string.format(
                "TRUE DEFECT %-14s initial=%.6e current=%.6e "
                    .. "relative=%.6e target=%.6e %s",
                componentNames[component], initialDefects[component],
                currentDefects[component], relative,
                flowAbsoluteTolerance,
                passed and "PASS" or "FAIL"
            ))
            if currentDefects[component] ~= currentDefects[component]
                or not passed then
                converged = false
            end
        end
        if converged then break end
    end
    assert(converged, label .. " restarted iterative/gmg"
        .. " failed the verified component tolerances")
    print("=== " .. label .. " OUTER solve converged ===")
    solution:set_consistent_storage_type()
    finalDefectField:set_consistent_storage_type()
    return solution, finalDefectField
end

coupledState, algebraicResidual = solveLinearSubproblem(
    approxSpace, coupledDisc, "Monolithic Stokes--Darcy"
)
coupledState:set_consistent_storage_type()
algebraicResidual:set_consistent_storage_type()

--------------------------------------------------------------------------------
-- Time discretization and initial values
--------------------------------------------------------------------------------
-- No time operator is assembled because the requested model is steady. The
-- zero algebraic initial guess is created inside solveLinearSubproblem. This
-- block converts the stored Darcy unknown to physical kinematic pressure for
-- all diagnostics and VTK fields; Stokes velocity and pressure are unchanged.

-- Convert only the stored Darcy unknown back to physical pressure. Using this
-- state for postprocessing also makes darcyEq:gradient() a physical gradient.
physicalState = coupledState:clone()
ScaleGF(physicalState, coupledState, {
    1.0, 1.0, 1.0, 1.0, darcyPressureScale
})
physicalState:set_consistent_storage_type()

-- A flow-only diagnostic uses the exact same refined mesh, equations, and
-- solver as the coupled model, then exits before any transport space or leak
-- operator is constructed.  This lets a new mesh be checked independently.
if flowOnly then
    local darcyMembraneVelocity = DarcyVelocityLinker()
    darcyMembraneVelocity:set_permeability(membranePermeability)
    darcyMembraneVelocity:set_viscosity(viscosity)
    darcyMembraneVelocity:set_density(1.0)
    darcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
    darcyMembraneVelocity:set_pressure_gradient(darcyMembraneEq:gradient())

    local darcyInterVelocity = DarcyVelocityLinker()
    darcyInterVelocity:set_permeability(interPermeability)
    darcyInterVelocity:set_viscosity(viscosity)
    darcyInterVelocity:set_density(1.0)
    darcyInterVelocity:set_gravity(ConstUserVector(0.0))
    darcyInterVelocity:set_pressure_gradient(darcyInterEq:gradient())

    local function vectorComponents(vectorData)
        local components = {}
        for component = 0, 2 do
            components[component + 1] = UserVectorEntryAdapter3d()
            components[component + 1]:set_vector(vectorData, component)
        end
        return components
    end

    local stokesOut = VTKOutput()
    stokesOut:clear_selection()
    stokesOut:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
    stokesOut:select_nodal("u", "u")
    stokesOut:select_nodal("v", "v")
    stokesOut:select_nodal("w", "w")
    stokesOut:select_nodal("p_stokes", "stokes_pressure_kinematic")
    stokesOut:print_subsets(
        vtk_file_name .. "_lumen", physicalState, "Lumen", 0, 0.0
    )

    local function writeDarcyOutput(suffix, subset, velocity)
        local components = vectorComponents(velocity)
        local output = VTKOutput()
        output:clear_selection()
        output:select_nodal(
            "p_darcy_scaled", "darcy_pressure_kinematic"
        )
        output:select_element(velocity, "darcy_velocity_mm_s")
        output:select_element(components[1], "u")
        output:select_element(components[2], "v")
        output:select_element(components[3], "w")
        output:print_subsets(
            vtk_file_name .. suffix, physicalState, subset, 0, 0.0
        )
    end
    writeDarcyOutput("_membrane", "Membrane", darcyMembraneVelocity)
    writeDarcyOutput("_inter", "Inter", darcyInterVelocity)

    local total_elapsed_time = GetClockS() - total_start_time
    print("FLOW-ONLY run complete; sodium transport was not constructed.")
    print(string.format(
        "Total flow-only runtime: %.6f mins", total_elapsed_time / 60.0
    ))
    return
end

--------------------------------------------------------------------------------
-- Sodium transport driven by the converged flow
--------------------------------------------------------------------------------
-- Sodium does not alter viscosity, permeability, or pressure, so this is an
-- exact one-way block coupling: solve Stokes--Darcy once, freeze its velocity,
-- then solve the linear advection--diffusion/leak system.  Explicit grid-
-- function data are essential here; they make velocity a known coefficient
-- instead of creating a bilinear concentration--flow term at the zero guess.

transportSpace = ApproximationSpace(dom)
transportSpace:add_fct(
    "na_lumen", "Lagrange", 1,
    "Lumen,Apical,Inlet,Outlet,TripleJunction"
)
transportSpace:add_fct(
    "na_membrane", "Lagrange", 1,
    "Membrane,Apical,Basolateral,MembraneWall,TripleJunction"
)
transportSpace:add_fct(
    "na_inter", "Lagrange", 1,
    "Inter,Basolateral,OuterWall,TripleJunction"
)
transportSpace:init_levels()
transportSpace:init_top_surface()

local solvedStokesVelocity = ExplicitGridFunctionVector(
    coupledState, "u,v,w"
)
local solvedScaledDarcyGradient = ExplicitGridFunctionGradient(
    coupledState, "p_darcy_scaled"
)

local transportDarcyMembraneVelocity = DarcyVelocityLinker()
transportDarcyMembraneVelocity:set_permeability(
    membranePermeability * darcyPressureScale
)
transportDarcyMembraneVelocity:set_viscosity(viscosity)
transportDarcyMembraneVelocity:set_density(1.0)
transportDarcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
transportDarcyMembraneVelocity:set_pressure_gradient(
    solvedScaledDarcyGradient
)

local transportDarcyInterVelocity = DarcyVelocityLinker()
transportDarcyInterVelocity:set_permeability(
    interPermeability * darcyPressureScale
)
transportDarcyInterVelocity:set_viscosity(viscosity)
transportDarcyInterVelocity:set_density(1.0)
transportDarcyInterVelocity:set_gravity(ConstUserVector(0.0))
transportDarcyInterVelocity:set_pressure_gradient(
    solvedScaledDarcyGradient
)

local naLumenEq = ConvectionDiffusionFV1("na_lumen", "Lumen")
naLumenEq:set_mass_scale(0.0)
naLumenEq:set_diffusion(naDiffusion)
naLumenEq:set_velocity(solvedStokesVelocity)
naLumenEq:set_upwind(FullUpwind())

local naMembraneEq = ConvectionDiffusionFV1("na_membrane", "Membrane")
naMembraneEq:set_mass_scale(0.0)
naMembraneEq:set_diffusion(naDiffusion)
naMembraneEq:set_velocity(transportDarcyMembraneVelocity)
naMembraneEq:set_upwind(FullUpwind())

local naInterEq = ConvectionDiffusionFV1("na_inter", "Inter")
naInterEq:set_mass_scale(0.0)
naInterEq:set_diffusion(naDiffusion)
naInterEq:set_velocity(transportDarcyInterVelocity)
naInterEq:set_upwind(FullUpwind())

local naBasolateralLeak = Leak({"na_membrane", "na_inter"})
naBasolateralLeak:set_permeability(naBasolateralPermeability)
local naBasolateralTransport = MembraneTransportFV1(
    "Basolateral", naBasolateralLeak
)
naBasolateralTransport:set_density_function(1.0)

local naBoundary = DirichletBoundary()
naBoundary:add(naInlet, "na_lumen", "Inlet")
naBoundary:add(naOuter, "na_inter", "OuterWall")
-- No Inter volume touches this endpoint support ring, so it carries no
-- physical control volume and is fixed to the prescribed outer concentration.
naBoundary:add(naOuter, "na_inter", "TripleJunction")

transportDisc = DomainDiscretization(transportSpace)
transportDisc:add(naLumenEq)
transportDisc:add(naMembraneEq)
transportDisc:add(naInterEq)
-- The apical Na membrane operator is deliberately absent for this diagnostic.
-- There is no Leak object and no MembraneTransportFV1 registered on Apical.
transportDisc:add(naBasolateralTransport)
transportDisc:add(naBoundary)
transportDisc:add(OneSideP1Constraints())

local function solveTransport()
    util.solver.defaults.approxSpace = transportSpace
    local matrix = AssembledLinearOperator(transportDisc)
    local rhs = GridFunction(transportSpace)
    local state = GridFunction(transportSpace)
    state:set(0.0)
    transportDisc:adjust_solution(state)
    transportDisc:assemble_linear(matrix, rhs)
    local emptyRows = RegularizeEmptyMatrixRows(matrix)
    if emptyRows > 0 then
        print(string.format(
            "Regularized %d empty shadow-surface transport rows.", emptyRows
        ))
    end
    local transportPreconditioner = ILUT(1.0e-8)
    transportPreconditioner:set_sort(true)
    -- The upwind advection and two interface leaks make this operator strongly
    -- nonsymmetric. A wider-restart GMRES is substantially more robust than
    -- BiCGStab. ILUT remains an iterative sparse preconditioner; no SuperLU or
    -- other direct solve is used. The larger Stokes--Darcy block above retains
    -- the geometric multigrid hierarchy requested for large simulations.
    local solver = GMRES(100)
    solver:set_preconditioner(transportPreconditioner)
    solver:set_convergence_check(
        ConvCheck(500, 1.0e-12, 1.0e-9, true)
    )
    solver:init(matrix, state)
    assert(solver:apply(state, rhs),
        "Sodium advection--diffusion/leak iterative solve failed")
    state:set_consistent_storage_type()

    local defect = rhs:clone()
    VecAssign(defect, rhs)
    matrix:apply_sub(defect, state)
    local factors = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    }
    local names = {"na_lumen", "na_membrane", "na_inter"}
    for component = 1, 3 do
        local part = defect:clone()
        ScaleGF(part, defect, factors[component])
        local norm = VecNorm(part)
        print(string.format(
            "TRUE TRANSPORT DEFECT %-11s abs=%.6e %s",
            names[component], norm,
            norm <= 1.0e-12 and "PASS" or "FAIL"
        ))
        assert(norm <= 1.0e-12,
            names[component] .. " transport residual is too large")
    end
    return state, defect
end

transportState, transportResidual = solveTransport()

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
-- DarcyVelocityLinker evaluates q_D=-(K/nu)grad(p_D) separately with membrane
-- and interstitial permeability. The three VTK datasets deliberately contain
-- physical pressure and velocity values, not the internally scaled unknown.
-- Integral diagnostics report pressure continuity, global lumen mass balance,
-- and the fraction leaving through Outlet versus Apical leakage.

darcyMembraneVelocity = DarcyVelocityLinker()
darcyMembraneVelocity:set_permeability(membranePermeability)
darcyMembraneVelocity:set_viscosity(viscosity)
darcyMembraneVelocity:set_density(1.0)
darcyMembraneVelocity:set_gravity(ConstUserVector(0.0))
darcyMembraneVelocity:set_pressure_gradient(darcyMembraneEq:gradient())

darcyInterVelocity = DarcyVelocityLinker()
darcyInterVelocity:set_permeability(interPermeability)
darcyInterVelocity:set_viscosity(viscosity)
darcyInterVelocity:set_density(1.0)
darcyInterVelocity:set_gravity(ConstUserVector(0.0))
darcyInterVelocity:set_pressure_gradient(darcyInterEq:gradient())

local pressureTraceJumpL2 = L2Error(
    physicalState, "p_stokes", physicalState, "p_darcy_scaled", 2,
    "Apical"
)
local stokesInletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Inlet", "Lumen"
)
local stokesOutletFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Outlet", "Lumen"
)
local stokesApicalFlow = IntegralNormalComponentOnManifold(
    stokesVelocity, physicalState, "Apical", "Lumen"
)
local darcyApicalFlow = IntegralNormalComponentOnManifold(
    darcyMembraneVelocity, physicalState, "Apical", "Membrane"
)
local inletPressureMean = Integral(
    physicalState, "p_stokes", "Inlet"
) / Integral(1.0, physicalState, "Inlet")
local outletPressureMean = Integral(
    physicalState, "p_stokes", "Outlet"
) / Integral(1.0, physicalState, "Outlet")
print(string.format(
    "Apical pressure-trace jump L2 (including constrained endpoint rings): %.6e",
    pressureTraceJumpL2
))
if outerWallNoFlux then
    print(string.format(
        "Mean Stokes pressure: Inlet=%.6e Pa, Outlet=%.6e Pa; OuterWall=no-flux",
        density * 1.0e-6 * inletPressureMean,
        density * 1.0e-6 * outletPressureMean
    ))
else
    print(string.format(
        "Mean Stokes pressure: Inlet=%.6e Pa, Outlet=%.6e Pa; OuterWall=%.6e Pa",
        density * 1.0e-6 * inletPressureMean,
        density * 1.0e-6 * outletPressureMean, outerPressurePa
    ))
end
print(string.format(
    "Lumen outward flows: Inlet=%.6e, Outlet=%.6e, Apical=%.6e, balance=%.6e mm^3/s",
    stokesInletFlow, stokesOutletFlow, stokesApicalFlow,
    stokesInletFlow + stokesOutletFlow + stokesApicalFlow
))
local imposedInflow = -stokesInletFlow
local conservativeLeakFlow = imposedInflow - stokesOutletFlow
print(string.format(
    "Conservative flow split: Outlet=%.3f%%, leakage=%.3f%% of inlet",
    100.0 * stokesOutletFlow / imposedInflow,
    100.0 * conservativeLeakFlow / imposedInflow
))
print(string.format(
    "Raw Apical surface integral=%.3f%% of inlet; quadrature/discrete balance offset=%.6e mm^3/s",
    100.0 * stokesApicalFlow / imposedInflow,
    stokesInletFlow + stokesOutletFlow + stokesApicalFlow
))
print(string.format(
    "Apical diagnostic: Stokes=%.6e, raw Darcy-gradient=%.6e, sum=%.6e mm^3/s",
    stokesApicalFlow, darcyApicalFlow,
    stokesApicalFlow + darcyApicalFlow
))

local function vectorComponents(vectorData)
    -- VTK receives both the vector and named scalar components so ParaView can
    -- display arrows or inspect u, v, and w independently.
    local components = {}
    for component = 0, 2 do
        components[component + 1] = UserVectorEntryAdapter3d()
        components[component + 1]:set_vector(vectorData, component)
    end
    return components
end
local darcyMembraneVelocityComponents = vectorComponents(darcyMembraneVelocity)
local darcyInterVelocityComponents = vectorComponents(darcyInterVelocity)

stokesOut = VTKOutput()
stokesOut:clear_selection()
stokesOut:select_nodal({"u", "v", "w"}, "stokes_velocity_mm_s")
stokesOut:select_nodal("u", "u")
stokesOut:select_nodal("v", "v")
stokesOut:select_nodal("w", "w")
stokesOut:select_nodal("p_stokes", "stokes_pressure_kinematic")
stokesOut:print_subsets(
    vtk_file_name .. "_lumen", physicalState, "Lumen", 0, 0.0
)

local function writeDarcyOutput(suffix, subset, velocity, components)
    -- Darcy velocity is element data because it is reconstructed from the P1
    -- pressure gradient, which is elementwise constant for affine elements.
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal("p_darcy_scaled", "darcy_pressure_kinematic")
    output:select_element(velocity, "darcy_velocity_mm_s")
    output:select_element(components[1], "u")
    output:select_element(components[2], "v")
    output:select_element(components[3], "w")
    output:print_subsets(
        vtk_file_name .. suffix, physicalState, subset, 0, 0.0
    )
end
writeDarcyOutput(
    "_membrane", "Membrane", darcyMembraneVelocity,
    darcyMembraneVelocityComponents
)
writeDarcyOutput(
    "_inter", "Inter", darcyInterVelocity,
    darcyInterVelocityComponents
)

-- A concentration-only PVTU combines all three compartments into one
-- ParaView dataset while retaining the compartment VTUs as its pieces.
local naPrefix = vtk_file_name .. "_Na"
local naOutputs = {
    {suffix = "_lumen", subset = "Lumen", field = "na_lumen"},
    {suffix = "_membrane", subset = "Membrane", field = "na_membrane"},
    {suffix = "_inter", subset = "Inter", field = "na_inter"}
}
for _, item in ipairs(naOutputs) do
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal(item.field, "Na")
    output:print_subsets(
        naPrefix .. item.suffix, transportState, item.subset, 0, 0.0
    )
end

if ProcRank() == 0 then
    local masterName = naPrefix .. "_t0000.pvtu"
    local master = assert(io.open(masterName, "w"),
        "Cannot create combined Na PVTU: " .. masterName)
    local baseName = assert(naPrefix:match("([^/]+)$"))
    master:write('<?xml version="1.0"?>\n')
    master:write('<VTKFile type="PUnstructuredGrid" version="0.1" ',
                 'byte_order="LittleEndian">\n')
    master:write('  <PUnstructuredGrid GhostLevel="0">\n')
    master:write('    <PPoints>\n')
    master:write('      <PDataArray type="Float32" NumberOfComponents="3"/>\n')
    master:write('    </PPoints>\n')
    master:write('    <PPointData>\n')
    master:write('      <PDataArray type="Float32" Name="Na" ',
                 'NumberOfComponents="1"/>\n')
    master:write('    </PPointData>\n')
    for _, item in ipairs(naOutputs) do
        if NumProcs() == 1 then
            master:write(string.format(
                '    <Piece Source="%s%s_t0000.vtu"/>\n',
                baseName, item.suffix
            ))
        else
            for rank = 0, NumProcs() - 1 do
                master:write(string.format(
                    '    <Piece Source="%s%s_p%d_t0000.vtu"/>\n',
                    baseName, item.suffix, rank
                ))
            end
        end
    end
    master:write('  </PUnstructuredGrid>\n')
    master:write('</VTKFile>\n')
    master:close()
    print("Combined ParaView Na output: " .. masterName)
end

-- Spatial map of the final assembled residual r=b-Ax. These values diagnose
-- where the discrete equations are least satisfied; they are not an estimate
-- of the continuum discretization error. In ParaView, use the Magnitude of
-- momentum_equation_residual or color by either pressure-equation residual.
local function writeResidualOutput(suffix, subset)
    local output = VTKOutput()
    output:clear_selection()
    output:select_nodal(
        {"u", "v", "w"}, "momentum_equation_residual"
    )
    output:select_nodal("u", "u_equation_residual")
    output:select_nodal("v", "v_equation_residual")
    output:select_nodal("w", "w_equation_residual")
    output:select_nodal(
        "p_stokes", "stokes_continuity_equation_residual"
    )
    output:select_nodal(
        "p_darcy_scaled", "darcy_equation_residual_scaled"
    )
    output:print_subsets(
        vtk_file_name .. suffix, algebraicResidual, subset, 0, 0.0
    )
end
writeResidualOutput("_residual_lumen", "Lumen")
writeResidualOutput("_residual_membrane", "Membrane")
writeResidualOutput("_residual_inter", "Inter")

local solve_elapsed_time = GetClockS() - solve_start_time
local total_elapsed_time = GetClockS() - total_start_time
print(string.format(
    "Monolithic Stokes--Darcy solve and postprocessing: %.6f mins",
    solve_elapsed_time / 60.0
))
print(string.format(
    "Total script runtime:                              %.6f mins",
    total_elapsed_time / 60.0
))
