-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/mesh_generation/archive/.*$"),
    "Expected this archived script below Model/mesh_generation/archive: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Generate a projector-backed curved lumen embedded in an interstitial box.
-- The centerline is a 20 mm quarter circle: inlet tangent +x, outlet tangent +y.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

local centerline = "curved_long_tube_centerline.swc"
local output = "runs/curved_long_tube.ugx"
local radius = 1.0
local anisotropy = 0.5
local boxPadding = 5.0
local bendRadius = 20.0 / (0.5 * math.pi)
local capTolerance = 0.05

local function file_contains(path, needle)
    local file = assert(io.open(path, "r"), "Cannot open: " .. path)
    local text = file:read("*a")
    file:close()
    return text:find(needle, 1, true) ~= nil
end

print("Generating curved long tube with ProMesh/neuro_collection")
print("  centerline: " .. centerline)
print("  output:     " .. output)
print("  radius:     " .. radius .. " mm")
print("  padding:    " .. boxPadding .. " mm")

-- This importer creates conforming Lumen and Inter tetrahedra inside a padded
-- OuterWall box and writes the NeuriteProjector/npSurfParams data needed for
-- geometry-following refinement.
import_neurites_with_inter_from_swc(
    centerline, output, anisotropy, 0, boxPadding
)

-- The importer initially labels the complete closed lumen surface Apical.
-- Split the two terminal caps with ProMesh. Cap vertices must share the
-- Inlet/Outlet subset so nodal Dirichlet data reaches the cap center and rim.
-- The rim edges remain Apical, which preserves projected wall-edge refinement.
local function label_terminal_caps(path)
    local mesh = Mesh()
    assert(LoadMesh(mesh, path), "ProMesh could not reload: " .. path)

    ClearSelection(mesh)
    SelectElementsByRangeX(
        mesh, -capTolerance, capTolerance, false, false, true, false
    )
    RestrictSelectionToSubset(mesh, 0) -- Apical
    AssignNewSubset(mesh, "Inlet", true, false, true, false)

    ClearSelection(mesh)
    SelectElementsByRangeY(
        mesh, bendRadius - capTolerance, bendRadius + capTolerance,
        false, false, true, false
    )
    RestrictSelectionToSubset(mesh, 0) -- remaining Apical faces
    AssignNewSubset(mesh, "Outlet", true, false, true, false)

    ClearSelection(mesh)
    AssignSubsetColors(mesh)
    assert(SaveMesh(mesh, path), "ProMesh could not save: " .. path)
end

-- Only the coarse mesh is written. The flow application loads this mesh and
-- constructs its projector-backed multigrid hierarchy in memory.
label_terminal_caps(output)

assert(file_contains(output, "<projector type=\"neurite\""),
    "Output is missing its serialized NeuriteProjector")
assert(file_contains(output, "<vertex_attachment name=\"npSurfParams\""),
    "Output is missing npSurfParams")
assert(file_contains(output, "<subset name=\"Lumen\""),
    "Output is missing the Lumen subset")
assert(file_contains(output, "<subset name=\"Inter\""),
    "Output is missing the Inter subset")
assert(file_contains(output, "<subset name=\"Apical\""),
    "Output is missing the lumen/inter interface subset")
assert(file_contains(output, "<subset name=\"OuterWall\""),
    "Output is missing the OuterWall subset")
assert(file_contains(output, "<subset name=\"Inlet\""),
    "Output is missing the Inlet subset")
assert(file_contains(output, "<subset name=\"Outlet\""),
    "Output is missing the Outlet subset")

print("Curved long-tube mesh generated successfully.")
print("The projected refinement hierarchy will be built in memory by the flow script.")
