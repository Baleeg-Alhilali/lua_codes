-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/.*$"),
    "Expected this script below Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- EXPERIMENTAL, NOT A VERIFIED SOLVER.
-- UG4's GlobalGridFunctionNumberData cannot locate every Apical quadrature
-- point on this projector mesh, so transfer of Stokes velocity to Darcy fails.
-- Retained as a reviewable partitioned-coupling prototype, not for results.
-- Steady Stokes--Darcy coupling on the projector-refined nephron grid.
-- Coordinates and velocities are in mm and s. Both pressures use the
-- kinematic unit mm^2/s^2; Pa = rho[kg/m^3] * 1e-6 * p_kinematic.
--
-- The Apical interface is solved by relaxed Dirichlet--Neumann iteration:
-- Darcy receives the Stokes normal velocity, while Stokes receives Darcy
-- pressure as normal traction and BJS tangential friction. All three terms
-- are updated until both fields converge. This is a partitioned solve.

ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
AssertPluginsLoaded({"NavierStokes", "ConvectionDiffusion", "SuperLU6"})

local startTime = GetClockS()
local dim = 3
local gridName = util.GetParam("-grid", "runs/third_trace_mm.ugx")
local centerlineName = util.GetParam(
    "-swc", "runs/third_trace_mm_outer_radius_1p5.swc")
local numRefs = math.floor(util.GetParamNumber("-numRefs", 0))
local nu = util.GetParamNumber("-visc", 1.0)          -- mm^2/s
local inflow = util.GetParamNumber("-inflow", 8.4)    -- mm/s
local density = util.GetParamNumber("-density", 1000.0) -- kg/m^3
local permeabilityM2 = util.GetParamNumber("-permeabilityM2", 1.0e-12)
local kMm2 = permeabilityM2 * 1.0e6
local alpha = util.GetParamNumber("-bjsAlpha", 0.2)
local outerPressurePa = util.GetParamNumber("-outerPressurePa", 0.0)
local relaxation = util.GetParamNumber("-couplingRelaxation", 0.2)
local tolerance = util.GetParamNumber("-couplingTolerance", 1.0e-5)
local maxIterations = math.floor(util.GetParamNumber("-maxCouplingIters", 30))
local outputPrefix = util.GetParam(
    "-vtk", "Results/third_trace_mm_stokes_darcy_bjs_ref" .. numRefs)

assert(numRefs >= 0 and nu > 0 and density > 0 and permeabilityM2 > 0)
assert(alpha > 0 and relaxation > 0 and relaxation <= 1)
assert(tolerance > 0 and maxIterations >= 1)

-- BJS law: (sigma_s n)_t = -(alpha*nu/sqrt(K)) u_s,t, in kinematic units.
local bjsFriction = alpha * nu / math.sqrt(kMm2) -- mm/s
local outerPressureKin = outerPressurePa / (density * 1.0e-6)
print(string.format(
    "K=%.6e mm^2, BJS friction=%.6e mm/s, OuterWall=%.6e Pa",
    kMm2, bjsFriction, outerPressurePa))

local function fileContains(path, needle)
    local f = assert(io.open(path, "r"), "Cannot open " .. path)
    for line in f:lines() do
        if line:find(needle, 1, true) then f:close(); return true end
    end
    f:close()
    return false
end

assert(fileContains(gridName, '<projector type="neurite"'))
assert(fileContains(gridName, 'name="npSurfParams"'))

local points = {}
for line in io.lines(centerlineName) do
    local a, b, c, d, x, y, z = line:match(
        "^%s*(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)")
    -- SWC columns are id, type, x, y, z, radius, parent.
    if a and tonumber(a) and tonumber(b) then
        points[#points + 1] = {tonumber(c), tonumber(d), tonumber(x)}
    end
end
assert(#points >= 2, "Missing SWC centerline")

-- Unit normal of the continuous projector surface, from Lumen to Membrane.
local function apicalNormal(x, y, z)
    local best = math.huge
    local rxBest, ryBest, rzBest = 0, 0, 0
    for i = 1, #points - 1 do
        local a, b = points[i], points[i + 1]
        local dx, dy, dz = b[1]-a[1], b[2]-a[2], b[3]-a[3]
        local len2 = dx*dx + dy*dy + dz*dz
        if len2 > 0 then
            local t = ((x-a[1])*dx + (y-a[2])*dy + (z-a[3])*dz) / len2
            t = math.max(0, math.min(1, t))
            local rx, ry, rz = x-a[1]-t*dx, y-a[2]-t*dy, z-a[3]-t*dz
            local r2 = rx*rx + ry*ry + rz*rz
            if r2 < best then
                best, rxBest, ryBest, rzBest = r2, rx, ry, rz
            end
        end
    end
    assert(best > 1.0e-24, "Apical normal undefined at centerline")
    local inv = 1.0 / math.sqrt(best)
    return rxBest*inv, ryBest*inv, rzBest*inv
end

InitUG(dim, AlgebraType("CPU", 1))
assert(NumProcs() == 1, "Partitioned boundary evaluation currently requires one MPI rank")
local preload = Domain()
preload:create_additional_subset_handler("projSH")
LoadDomain(preload, gridName) -- registers npSurfParams before the actual load
preload = nil
collectgarbage("collect")
local dom = Domain()
dom:create_additional_subset_handler("projSH")
LoadDomain(dom, gridName)
assert(util.CheckSubsets(dom, {
    "Lumen", "Membrane", "Inter", "Apical", "Basolateral",
    "OuterWall", "Inlet", "Outlet"}))
balancer.RefineAndRebalanceDomain(dom, numRefs)
print(dom:domain_info():to_string())

local stokesSpace = ApproximationSpace(dom)
for _, name in ipairs({"u", "v", "w", "p"}) do
    stokesSpace:add_fct(name, "Lagrange", 1, "Lumen,Apical,Inlet,Outlet")
end
stokesSpace:init_levels()
stokesSpace:init_top_surface()

local darcySpace = ApproximationSpace(dom)
darcySpace:add_fct(
    "q", "Lagrange", 1,
    "Membrane,Inter,Apical,Basolateral,OuterWall,Inlet,Outlet")
darcySpace:init_levels()
darcySpace:init_top_surface()

local stokesState = GridFunction(stokesSpace)
local darcyState = GridFunction(darcySpace)
stokesState:set(0.0)
darcyState:set(0.0)
local qData = GlobalGridFunctionNumberData(darcyState, "q")
local uxData = GlobalGridFunctionNumberData(stokesState, "u")
local uyData = GlobalGridFunctionNumberData(stokesState, "v")
local uzData = GlobalGridFunctionNumberData(stokesState, "w")

-- Step inside the corresponding volume to avoid point-location ambiguity at
-- the projector surface. 0.01 um is much smaller than the 1.5 um lumen radius.
local sampleOffset = 0.0
-- UG4 probes Lua callbacks with (0,0,0) when registering a boundary.
local callbacksReady = false
local hasStokes = false
local hasDarcy = false
local function sample(data, x, y, z, nx, ny, nz, side)
    return data:evaluate_global({
        x+side*sampleOffset*nx,
        y+side*sampleOffset*ny,
        z+side*sampleOffset*nz})
end
local function oldVelocity(x, y, z, nx, ny, nz)
    if not hasStokes then return 0.0, 0.0, 0.0 end
    return sample(uxData,x,y,z,nx,ny,nz,-1),
           sample(uyData,x,y,z,nx,ny,nz,-1),
           sample(uzData,x,y,z,nx,ny,nz,-1)
end
local function oldDarcyPressure(x, y, z, nx, ny, nz)
    if not hasDarcy then return 0.0 end
    return sample(qData,x,y,z,nx,ny,nz,1)
end

function InletVelocity(x, y, z, t)
    return -0.771159*inflow, -0.349961*inflow, 0.524920*inflow
end

-- NeumannBoundaryFV1 takes outward diffusive flux, which is the negative
-- of prescribed physical traction in the Stokes momentum weak form.
-- Its vector callback is dotted with the discrete face normal, giving q*n_i.
function DarcyPressureFluxU(x, y, z, t, si)
    if not callbacksReady then return 0, 0, 0 end
    local nx, ny, nz = apicalNormal(x,y,z)
    return oldDarcyPressure(x,y,z,nx,ny,nz), 0, 0
end
function DarcyPressureFluxV(x, y, z, t, si)
    if not callbacksReady then return 0, 0, 0 end
    local nx, ny, nz = apicalNormal(x,y,z)
    return 0, oldDarcyPressure(x,y,z,nx,ny,nz), 0
end
function DarcyPressureFluxW(x, y, z, t, si)
    if not callbacksReady then return 0, 0, 0 end
    local nx, ny, nz = apicalNormal(x,y,z)
    return 0, 0, oldDarcyPressure(x,y,z,nx,ny,nz)
end

local function bjsFlux(component, x, y, z)
    local nx, ny, nz = apicalNormal(x,y,z)
    local ux, uy, uz = oldVelocity(x,y,z,nx,ny,nz)
    local un = ux*nx + uy*ny + uz*nz
    local u = ({ux,uy,uz})[component]
    local n = ({nx,ny,nz})[component]
    return bjsFriction * (u - un*n)
end
function BJSFluxU(x,y,z,t,si)
    if not callbacksReady then return 0 end
    return bjsFlux(1,x,y,z)
end
function BJSFluxV(x,y,z,t,si)
    if not callbacksReady then return 0 end
    return bjsFlux(2,x,y,z)
end
function BJSFluxW(x,y,z,t,si)
    if not callbacksReady then return 0 end
    return bjsFlux(3,x,y,z)
end

-- The FV1 vector Neumann operator uses the actual Membrane face normal.
-- The same velocity is sent across Apical, so normal volume flux is matched.
function StokesVelocityForDarcy(x,y,z,t,si)
    if not callbacksReady then return 0, 0, 0 end
    local nx, ny, nz = apicalNormal(x,y,z)
    return oldVelocity(x,y,z,nx,ny,nz)
end

local ns = NavierStokesFV1({"u","v","w","p"}, {"Lumen"})
ns:set_stokes(true)
ns:set_laplace(false) -- symmetric-gradient stress for BJS traction
ns:set_kinematic_viscosity(nu)
ns:set_exact_jacobian(true)
ns:set_stabilization("fields", "cor")
local inlet = NavierStokesInflow(ns)
inlet:add("InletVelocity", "Inlet")
local openBoundary = NavierStokesNoNormalStressOutflow(ns)
openBoundary:add("Outlet,Apical")

local stokesDisc = DomainDiscretization(stokesSpace)
stokesDisc:add(ns)
stokesDisc:add(inlet)
stokesDisc:add(openBoundary)
for _, row in ipairs({
    {"u", "DarcyPressureFluxU", "BJSFluxU"},
    {"v", "DarcyPressureFluxV", "BJSFluxV"},
    {"w", "DarcyPressureFluxW", "BJSFluxW"}}) do
    local bc = NeumannBoundary(row[1], "fv1")
    bc:add(row[2], "Apical", "Lumen")
    bc:add(row[3], "Apical", "Lumen")
    stokesDisc:add(bc)
end

local membraneEq = ConvectionDiffusionFV1("q", "Membrane")
membraneEq:set_mass_scale(0.0)
membraneEq:set_diffusion(kMm2 / nu)
local interEq = ConvectionDiffusionFV1("q", "Inter")
interEq:set_mass_scale(0.0)
interEq:set_diffusion(kMm2 / nu)
local apicalFlux = NeumannBoundary("q", "fv1")
apicalFlux:add("StokesVelocityForDarcy", "Apical", "Membrane")
local outerPressure = DirichletBoundary()
outerPressure:add(outerPressureKin, "q", "OuterWall")
local darcyDisc = DomainDiscretization(darcySpace)
darcyDisc:add(membraneEq)
darcyDisc:add(interEq)
darcyDisc:add(apicalFlux)
darcyDisc:add(outerPressure)
callbacksReady = true

-- Both subproblems are linear for frozen interface data. Under-relaxation
-- damps the pressure/velocity feedback on this thin, high-resistance layer.
local stokesOld = GridFunction(stokesSpace)
local darcyOld = GridFunction(darcySpace)
local stokesDelta = GridFunction(stokesSpace)
local darcyDelta = GridFunction(darcySpace)
local converged = false
for iter = 1, maxIterations do
    VecScaleAdd2(stokesOld, 1.0, stokesState, 0.0, stokesState)
    util.solver.defaults.approxSpace = stokesSpace
    local _, newStokes = util.solver.SolveLinearProblem(stokesDisc, "superlu")
    VecScaleAdd2(stokesState, relaxation, newStokes,
                 1.0-relaxation, stokesOld)
    hasStokes = true

    VecScaleAdd2(darcyOld, 1.0, darcyState, 0.0, darcyState)
    util.solver.defaults.approxSpace = darcySpace
    local _, newDarcy = util.solver.SolveLinearProblem(darcyDisc, "superlu")
    VecScaleAdd2(darcyState, relaxation, newDarcy,
                 1.0-relaxation, darcyOld)
    hasDarcy = true

    VecScaleAdd2(stokesDelta, 1.0, stokesState, -1.0, stokesOld)
    VecScaleAdd2(darcyDelta, 1.0, darcyState, -1.0, darcyOld)
    local relU = VecNorm(stokesDelta) / math.max(1.0, VecNorm(stokesState))
    local relQ = VecNorm(darcyDelta) / math.max(1.0, VecNorm(darcyState))
    print(string.format(
        "Coupling iteration %d: relative velocity %.4e, pressure %.4e",
        iter, relU, relQ))
    if math.max(relU, relQ) < tolerance then
        converged = true
        break
    end
end
assert(converged, "Stokes--Darcy interface iteration did not converge")

local membraneVelocity = DarcyVelocityLinker()
membraneVelocity:set_permeability(kMm2)
membraneVelocity:set_viscosity(nu)
membraneVelocity:set_density(1.0)
membraneVelocity:set_gravity(ConstUserVector(0.0))
membraneVelocity:set_pressure_gradient(membraneEq:gradient())
local interVelocity = DarcyVelocityLinker()
interVelocity:set_permeability(kMm2)
interVelocity:set_viscosity(nu)
interVelocity:set_density(1.0)
interVelocity:set_gravity(ConstUserVector(0.0))
interVelocity:set_pressure_gradient(interEq:gradient())

local lumenOut = VTKOutput()
lumenOut:clear_selection()
lumenOut:select_nodal({"u","v","w"}, "stokes_velocity_mm_s")
lumenOut:select_nodal("p", "stokes_pressure_kinematic")
lumenOut:print_subsets(outputPrefix .. "_lumen", stokesState, "Lumen", 0, 0.0)
local membraneOut = VTKOutput()
membraneOut:clear_selection()
membraneOut:select_nodal("q", "darcy_pressure_kinematic")
membraneOut:select_element(membraneVelocity, "darcy_velocity_mm_s")
membraneOut:print_subsets(outputPrefix .. "_membrane", darcyState, "Membrane", 0, 0.0)
local interOut = VTKOutput()
interOut:clear_selection()
interOut:select_nodal("q", "darcy_pressure_kinematic")
interOut:select_element(interVelocity, "darcy_velocity_mm_s")
interOut:print_subsets(outputPrefix .. "_inter", darcyState, "Inter", 0, 0.0)

print(string.format("Coupled solve completed in %.3f min",
    (GetClockS()-startTime)/60.0))
