-- Resolve every relative input and output path from the Model directory.
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

total_start_time = GetClockS()

--------------------------------------------------------------------------------
-- Geometry and time interval
--------------------------------------------------------------------------------

dim = 2
gridS = "2D_trans2"
gridName = "mesh_generation/meshes/legacy_transport_meshes/" .. gridS .. ".ugx"


numRefs = util.GetParamNumber("-numRefs", 4)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)
vtk_file_name = "Results_final/"..gridS .. numRefs
endTime = util.GetParamNumber("-endTime", 15.0) -- [s]
dt = util.GetParamNumber("-dt", 0.1)           -- [s]

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- The mesh coordinates are in micrometres, so diffusion is in um^2/s.
D_na = 1.33e-1

-- Leak permeability must use the same length unit as the mesh.
m_to_um = 1.0e6
P_apical = 1.0e-6 * m_to_um       -- [um/s], Lumen -> Membrane
P_basolateral = 5.0e-7 * m_to_um  -- [um/s], Membrane -> Inter

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

InitUG(dim, AlgebraType("CPU", 1))
print("MPI processes: " .. NumProcs())

--------------------------------------------------------------------------------
-- Domain and approximation space
--------------------------------------------------------------------------------

dom = util.CreateDomain(
    gridName,
    numPreRefs,
    {
        "Lumen",
        "Membrane",
        "Basolateral",
        "Apical",
        "Inter",
        "Toplum",
        "Topmem",
        "Topint",
        "LeftBnd",
        "RightBnd"
    },
    true
)

balancer.RefineAndRebalanceDomain(dom, numRefs - numPreRefs)

print("Domain info:")
print(dom:domain_info():to_string())

-- Separate concentration functions are required. Using a single function over
-- all three volumes would make the concentration continuous across the
-- interfaces and bypass the membrane resistance.
approxSpace = ApproximationSpace(dom)

approxSpace:add_fct(
    "c_lumen",
    "Lagrange",
    1,
    "Lumen,Apical,LeftBnd,Toplum"
)

approxSpace:add_fct(
    "c_membrane",
    "Lagrange",
    1,
    "Membrane,Apical,Basolateral,Topmem"
)

approxSpace:add_fct(
    "c_inter",
    "Lagrange",
    1,
    "Inter,Basolateral,RightBnd,Topint"
)

approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "x")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

-- Apical is now a membrane-transport interface. Therefore, the prescribed
-- lumen concentration is applied at the inlet instead.
dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c_lumen", "LeftBnd")
dirichletBND:add(0.0, "c_inter", "RightBnd")

-- Diffusion inside the lumen.
lumenDiffusion = ConvectionDiffusionFV1("c_lumen", "Lumen")
lumenDiffusion:set_diffusion(D_na)

-- Diffusion inside the membrane volume.
membraneDiffusion = ConvectionDiffusionFV1(
    "c_membrane",
    "Membrane"
)
membraneDiffusion:set_diffusion(D_na)

-- Diffusion inside the interstitial volume.
interDiffusion = ConvectionDiffusionFV1("c_inter", "Inter")
interDiffusion:set_diffusion(D_na)

--------------------------------------------------------------------------------
-- Apical membrane transport: Lumen -> Membrane
--------------------------------------------------------------------------------

-- Leak calculates:
--
--     flux = permeability * (source - target)
--
-- The order of the concentration functions therefore sets the positive
-- transport direction.
apicalLeak = Leak({
    "c_lumen",
    "c_membrane"
})

apicalLeak:set_permeability(P_apical)

apicalTransport = MembraneTransportFV1(
    "Apical",
    apicalLeak
)

apicalTransport:set_density_function(0.4)

--------------------------------------------------------------------------------
-- Basolateral membrane transport: Membrane -> Inter
--------------------------------------------------------------------------------

basolateralLeak = Leak({
    "c_membrane",
    "c_inter"
})

basolateralLeak:set_permeability(P_basolateral)

basolateralTransport = MembraneTransportFV1(
    "Basolateral",
    basolateralLeak
)

basolateralTransport:set_density_function(1.0)

--------------------------------------------------------------------------------
-- Domain discretization
--------------------------------------------------------------------------------

domainDisc = DomainDiscretization(approxSpace)

domainDisc:add(lumenDiffusion)
domainDisc:add(membraneDiffusion)
domainDisc:add(interDiffusion)

domainDisc:add(apicalTransport)
domainDisc:add(basolateralTransport)

domainDisc:add(dirichletBND)

--------------------------------------------------------------------------------
-- MPI-friendly iterative solver
--------------------------------------------------------------------------------

util.solver.defaults.approxSpace = approxSpace

solverDesc = {
    type = "newton", -- alternative for linearized runs: a linear solver

    convCheck = {
        type = "standard",
        iterations = 5,
        absolute = 1e-14,
        reduction = 1e-8,
        verbose = true
    },

    linSolver = {
        type = "bicgstab", -- alternative: "gmres"

        precond = {
            type = "gmg", -- alternative: "ilu" for a single-level solve

            smoother = {
                type = "gs",
                overlap = true
            },

            cycle = "V",
            preSmooth = 4,
            postSmooth = 4,
            rap = true,
            baseLevel = 0,

            baseSolver = {
                type = "bicgstab",

                precond = {
                    type = "ilu",
                    overlap = true
                },

                convCheck = {
                    type = "standard",
                    iterations = 200,
                    absolute = 1e-14,
                    reduction = 0.5e-8,
                    verbose = false
                }
            }
        },

        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-14,
            reduction = 0.5e-8,
            verbose = true
        }
    }
}

solver = util.solver.CreateSolver(solverDesc)

u = GridFunction(approxSpace)

Interpolate(0.0, u, "c_lumen")
Interpolate(0.0, u, "c_membrane")
Interpolate(0.0, u, "c_inter")
--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
vtkLumen = VTKOutput()
vtkLumen:select("c_lumen", "Na")

vtkMembrane = VTKOutput()
vtkMembrane:select("c_membrane", "Na")

vtkInter = VTKOutput()
vtkInter:select("c_inter", "Na")

-- Each compartment is written as a VTU piece because its concentration field
-- is subset-restricted. One PVTU master combines all pieces for convenient
-- loading as a single dataset in ParaView.
local outputTimes = {}
local lastOutputStep = 0

local function writeCombinedPVTU(step, time)
    if ProcRank() ~= 0 then return end

    local fileName = string.format("%s_t%04d.pvtu", vtk_file_name, step)
    local file = assert(io.open(fileName, "w"),
        "Cannot create combined PVTU file: " .. fileName)
    local baseName = assert(vtk_file_name:match("([^/]+)$"))
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
    file:write('      <PDataArray type="Float32" Name="Na" ',
               'NumberOfComponents="1"/>\n')
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

function writeVTK(uOut, step, time)
    vtkLumen:print_subsets(
        vtk_file_name .. "_lumen",
        uOut,
        "Lumen",
        step,
        time
    )

    vtkMembrane:print_subsets(
        vtk_file_name .. "_membrane",
        uOut,
        "Membrane",
        step,
        time
    )

    vtkInter:print_subsets(
        vtk_file_name .. "_inter",
        uOut,
        "Inter",
        step,
        time
    )

    writeCombinedPVTU(step, time)
end

-- Write the initialized solution at t = 0.
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
timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
timeIntegrator:attach_observer(vtkObserver)

--solver:init(AssembledOperator(domainDisc))
--solver:prepare(u)

solve_start_time = GetClockS()
timeIntegrator:apply(
    u,
    endTime,
    u,
    0.0
)

-- Create one time-series index whose entries are the combined PVTU masters.
if ProcRank() == 0 then
    local pvdFileName = vtk_file_name .. ".pvd"
    local pvd = assert(io.open(pvdFileName, "w"),
        "Cannot create combined PVD file: " .. pvdFileName)
    local baseName = assert(vtk_file_name:match("([^/]+)$"))

    pvd:write('<?xml version="1.0"?>\n')
    pvd:write('<VTKFile type="Collection" version="0.1" ',
              'byte_order="LittleEndian">\n')
    pvd:write('  <Collection>\n')
    for step = 0, lastOutputStep do
        if outputTimes[step] ~= nil then
            pvd:write(string.format(
                '    <DataSet timestep="%.17g" group="" part="0" ' ..
                'file="%s_t%04d.pvtu"/>\n',
                outputTimes[step], baseName, step))
        end
    end
    pvd:write('  </Collection>\n')
    pvd:write('</VTKFile>\n')
    pvd:close()
end

solve_elapsed_time = GetClockS() - solve_start_time
total_elapsed_time = GetClockS() - total_start_time

print(
    string.format(
        "Time integration runtime: %.6f mins",
        solve_elapsed_time / 60.0
    )
)

print(
    string.format(
        "Total script runtime:      %.6f mins",
        total_elapsed_time / 60.0
    )
)
