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

dim = 3
gridS = "3rd_trace_third_boxref4"
gridName = "ProMeshFiles/runs/" .. gridS .. ".ugx"

numRefs = util.GetParamNumber("-numRefs", 2)
numPreRefs = util.GetParamNumber("-numPreRefs", 0)

endTime = util.GetParamNumber("-endTime", 0.020) -- [s]
dt = util.GetParamNumber("-dt", 0.001)           -- [s]

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

-- The mesh coordinates are in micrometres, so diffusion is in um^2/s.
D_na = 1.33e3

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
        "OuterWall",
        "Inlet",
        "Outlet"
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
    "Lumen"
)

approxSpace:add_fct(
    "c_membrane",
    "Lagrange",
    1,
    "Membrane"
)

approxSpace:add_fct(
    "c_inter",
    "Lagrange",
    1,
    "Inter"
)

approxSpace:init_levels()
approxSpace:init_top_surface()

print("Approximation space:")
approxSpace:print_statistic()

OrderLex(approxSpace, "xy")

--------------------------------------------------------------------------------
-- Spatial discretization
--------------------------------------------------------------------------------

-- Apical is now a membrane-transport interface. Therefore, the prescribed
-- lumen concentration is applied at the inlet instead.
dirichletBND = DirichletBoundary()
dirichletBND:add(1.0, "c_lumen", "Inlet")

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

apicalTransport:set_density_function(1.0)

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
    type = "newton",

    convCheck = {
        type = "standard",
        iterations = 5,
        absolute = 1e-8,
        reduction = 1e-8,
        verbose = true
    },

    linSolver = {
        type = "bicgstab",

        precond = {
            type = "gmg",

            smoother = {
                type = "gs",
                overlap = true
            },

            cycle = "V",
            preSmooth = 4,
            postSmooth = 4,
            rap = true,
            baseLevel = numPreRefs,

            baseSolver = {
                type = "bicgstab",

                precond = {
                    type = "ilu",
                    overlap = true
                },

                convCheck = {
                    type = "standard",
                    iterations = 200,
                    absolute = 1e-10,
                    reduction = 0.5e-8,
                    verbose = false
                }
            }
        },

        convCheck = {
            type = "standard",
            iterations = 100,
            absolute = 1e-8,
            reduction = 0.5e-8,
            verbose = true
        }
    }
}

solver = util.solver.CreateSolver(solverDesc)
--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------
--out = VTKOutput()
--out:clear_selection()
--out:select_nodal("c_lumen", "c_lumen")
--out:select_nodal("c_membrane", "c_membrane")
--out:select_nodal("c_inter", "c_inter")
vtkLumen = VTKOutput()
vtkLumen:select("c_na_lumen", "Na")
vtkLumen:print_subsets(
    "vtk/Na_lumen",
    u,
    "Lumen",
    step,
    time
)

vtkMembrane = VTKOutput()
vtkMembrane:select("c_na_membrane", "Na")
vtkMembrane:print_subsets(
    "vtk/Na_membrane",
    u,
    "Membrane",
    step,
    time
)

vtkInter = VTKOutput()
vtkInter:select("c_na_inter", "Na")
vtkInter:print_subsets(
    "vtk/Na_inter",
    u,
    "Inter",
    step,
    time
)
--------------------------------------------------------------------------------
-- Time stepping
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc, 1.0)

timeIntegrator = SimpleTimeIntegrator(timeDisc)
timeIntegrator:set_solver(solver)
timeIntegrator:set_time_step(dt)
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name.."1", vtkLumen))
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name.."2", vtkMembrane))
timeIntegrator:attach_observer(VTKOutputObserver(vtk_file_name.."3", vtkInter))



u = GridFunction(approxSpace)

Interpolate(0.0, u, "c_lumen")
Interpolate(0.0, u, "c_membrane")
Interpolate(0.0, u, "c_inter")

solver:init(AssembledOperator(domainDisc))
solver:prepare(u)

solve_start_time = GetClockS()

timeIntegrator:apply(u,endTime,u,0.0)

vtk:write_time_pvd(out, u)

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
