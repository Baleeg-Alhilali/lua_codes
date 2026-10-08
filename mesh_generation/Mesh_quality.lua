-- Resolve CLI paths from the mesh-generation directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local GENERATION_DIR = assert(__model_script:match("^(.*)/[^/]+$"),
    "Cannot determine mesh_generation directory: " .. __model_script)
local MODEL_ROOT = assert(GENERATION_DIR:match("^(.*)/mesh_generation$"),
    "Expected this script inside Model/mesh_generation: " .. __model_script)
ChangeDirectory(GENERATION_DIR)

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh"})

-- ============================================================================
-- Command-line parameters
-- ============================================================================

local inputFile = util.GetParam("-grid", "")
local histogramBins = math.floor(util.GetParamNumber("-bins", 20))
local sliverThreshold = util.GetParamNumber("-sliverThreshold", 0.01)
local anisotropyThreshold =
    util.GetParamNumber("-anisotropyThreshold", 0.1)

local saveSelections = util.HasParamOption("-saveSelections")

if inputFile == "" then
    error(
        "Missing input grid.\n" ..
        "Example:\n" ..
        "ugshell -ex Mesh_quality.lua " ..
        "-grid mesh.ugx -saveSelections"
    )
end

if histogramBins < 1 then
    error("-bins must be at least 1.")
end

if sliverThreshold <= 0.0 then
    error("-sliverThreshold must be greater than zero.")
end

if anisotropyThreshold <= 0.0 then
    error("-anisotropyThreshold must be greater than zero.")
end

-- ============================================================================
-- Helper functions
-- ============================================================================

local function separator(title)
    print("")
    print("============================================================")
    print(title)
    print("============================================================")
end

local function stripExtension(path)
    return string.gsub(path, "%.[^%.%/]+$", "")
end

local function countAllVolumes(obj)
    local iterator = obj:volume_begin()
    local iteratorEnd = obj:volume_end()
    local count = 0

    while iterator:unequal(iteratorEnd) do
        count = count + 1
        iterator:advance()
    end

    return count
end

local function countSelectedVolumes(obj)
    local iterator = obj:volume_selection_begin()
    local iteratorEnd = obj:volume_selection_end()
    local count = 0

    while iterator:unequal(iteratorEnd) do
        count = count + 1
        iterator:advance()
    end

    return count
end

local function percentage(count, total)
    if total == 0 then
        return 0.0
    end

    return 100.0 * count / total
end

local function printCount(label, count, total)
    print(
        string.format(
            "%s: %d of %d (%.6f%%)",
            label,
            count,
            total,
            percentage(count, total)
        )
    )
end

local function countVolumeType(
    obj,
    selectHexahedra,
    selectOctahedra,
    selectPrisms,
    selectPyramids,
    selectTetrahedra
)
    ClearSelection(obj)

    SelectVolumesByType(
        obj,
        selectHexahedra,
        selectOctahedra,
        selectPrisms,
        selectPyramids,
        selectTetrahedra
    )

    local count = countSelectedVolumes(obj)
    ClearSelection(obj)

    return count
end

-- ============================================================================
-- Load mesh
-- ============================================================================

separator("LOADING GRID")

print("Input file: " .. inputFile)

local mesh = Mesh()

if not LoadMesh(mesh, inputFile) then
    error("Could not load grid: " .. inputFile)
end

local outputPrefix = stripExtension(inputFile)
local totalVolumes = countAllVolumes(mesh)

print("Grid loaded successfully.")
print("Total volumes: " .. totalVolumes)

-- ============================================================================
-- Determine volume composition
-- ============================================================================

separator("VOLUME COMPOSITION")

local tetrahedronCount = countVolumeType(
    mesh,
    false,  -- hexahedra
    false,  -- octahedra
    false,  -- prisms
    false,  -- pyramids
    true    -- tetrahedra
)

local pyramidCount = countVolumeType(
    mesh,
    false,
    false,
    false,
    true,
    false
)

local prismCount = countVolumeType(
    mesh,
    false,
    false,
    true,
    false,
    false
)

local octahedronCount = countVolumeType(
    mesh,
    false,
    true,
    false,
    false,
    false
)

local hexahedronCount = countVolumeType(
    mesh,
    true,
    false,
    false,
    false,
    false
)

local recognizedVolumes =
    tetrahedronCount +
    pyramidCount +
    prismCount +
    octahedronCount +
    hexahedronCount

local otherVolumeCount = totalVolumes - recognizedVolumes

print("The grid contains the following volume-element types:")
print("")
print("  Tetrahedra:  " .. tetrahedronCount)
print("  Pyramids:    " .. pyramidCount)
print("  Prisms:      " .. prismCount)
print("  Octahedra:   " .. octahedronCount)
print("  Hexahedra:   " .. hexahedronCount)
print("  Other:       " .. otherVolumeCount)
print("  -------------------------")
print("  All volumes: " .. totalVolumes)

print("")
print([[
A tetrahedron has 4 vertices.
A pyramid has 5 vertices.
A prism has 6 vertices.
A hexahedron has 8 vertices.

The ProMesh aspect-ratio function used below supports tetrahedra only.
]])

-- ============================================================================
-- 1. Overall tetrahedron quality
-- ============================================================================

separator("1. OVERALL TETRAHEDRON QUALITY")

print([[
The ProMesh tetrahedron aspect ratio is normalized to the range 0 to 1:

    1.0       = an ideal regular tetrahedron
    Near 1.0  = a well-shaped tetrahedron
    Near 0.0  = a flat, stretched, sliver-like, or degenerate tetrahedron

The printed statistics mean:

    min  = quality of the worst tetrahedron
    max  = quality of the best tetrahedron
    mean = average quality of all tetrahedra
    sd   = standard deviation

A large standard deviation means that tetrahedron quality varies
considerably across the mesh.

Only tetrahedra are selected for this calculation because the ProMesh
aspect-ratio function does not support other volume types.
]])

ClearSelection(mesh)

SelectVolumesByType(
    mesh,
    false,  -- hexahedra
    false,  -- octahedra
    false,  -- prisms
    false,  -- pyramids
    true    -- tetrahedra
)

local selectedTetrahedronCount = countSelectedVolumes(mesh)

print("Selected tetrahedra: " .. selectedTetrahedronCount)
print("Total grid volumes:  " .. totalVolumes)

if selectedTetrahedronCount > 0 then
    print("")
    print("Tetrahedron aspect-ratio statistics:")
    print("")

    PrintVolumeAspectRatios(mesh)

    print("")
    print(
        "Aspect-ratio histogram with " ..
        histogramBins ..
        " bins:"
    )

    print([[
Each histogram row gives:

    aspect-ratio interval | number of tetrahedra in that interval

Bins near 0 contain poor-quality tetrahedra.
Bins near 1 contain better-quality tetrahedra.
]])

    PrintVolumeAspectRatioHistogram(mesh, histogramBins)
else
    print("")
    print("No tetrahedra were found.")
    print("The aspect-ratio calculation was skipped.")
end

ClearSelection(mesh)

-- ============================================================================
-- 2. Sliver tetrahedra
-- ============================================================================

separator("2. VERY FLAT / SLIVER TETRAHEDRA")

print(
    string.format(
        [[
The sliver-selection threshold is %.8g.

A sliver is a tetrahedron with very little height or volume relative to
its edge lengths. A sliver can have edges of apparently reasonable
length while enclosing almost no volume.

Slivers can cause:

    * poorly conditioned equation systems
    * inaccurate numerical gradients
    * unstable or slow solver convergence
    * unreliable finite-element or finite-volume results

The number reported below is the number of tetrahedra that ProMesh
classified as slivers.

A lower threshold selects only more severe slivers.
]],
        sliverThreshold
    )
)

ClearSelection(mesh)

local sliverCount = SelectSlivers(mesh, sliverThreshold)

printCount(
    "Selected sliver tetrahedra",
    sliverCount,
    tetrahedronCount
)

if sliverCount == 0 then
    print("Result: no sliver tetrahedra were detected.")
else
    print("Result: inspect the selected sliver tetrahedra.")
end

if saveSelections then
    local sliverFile = outputPrefix .. "_selected_slivers.ugx"

    if SaveMesh(mesh, sliverFile) then
        print("Saved sliver selection to: " .. sliverFile)
    else
        print("WARNING: could not save: " .. sliverFile)
    end
end

ClearSelection(mesh)

-- ============================================================================
-- 3. Anisotropic volumes
-- ============================================================================

separator("3. EXTREMELY ANISOTROPIC VOLUMES")

print(
    string.format(
        [[
The anisotropy threshold is %.8g.

This test compares short and long edge scales in each volume.

    Ratio near 1 = edge lengths are relatively similar
    Ratio near 0 = the volume is strongly stretched

With a threshold of %.8g, ProMesh selects volumes whose edge ratio is
below the threshold.

An anisotropic volume is not automatically invalid. Anisotropy may be
intentional when the mesh follows a long thin structure. However,
extreme anisotropy should be checked because it can affect numerical
accuracy and solver conditioning.
]],
        anisotropyThreshold,
        anisotropyThreshold
    )
)

ClearSelection(mesh)
SelectAnisotropicVolumes(mesh, anisotropyThreshold)

local anisotropicCount = countSelectedVolumes(mesh)

printCount(
    "Selected anisotropic volumes",
    anisotropicCount,
    totalVolumes
)

if anisotropicCount == 0 then
    print("Result: no volumes exceeded the anisotropy criterion.")
else
    print("Result: inspect whether the selected anisotropy is intentional.")
end

if saveSelections then
    local anisotropicFile =
        outputPrefix .. "_selected_anisotropic.ugx"

    if SaveMesh(mesh, anisotropicFile) then
        print("Saved anisotropic selection to: " .. anisotropicFile)
    else
        print("WARNING: could not save: " .. anisotropicFile)
    end
end

ClearSelection(mesh)

-- ============================================================================
-- 4. Unorientable/pathological volumes
-- ============================================================================

separator("4. UNORIENTABLE / PATHOLOGICAL VOLUMES")

print([[
An orientable volume has a geometrically meaningful signed orientation.

An unorientable volume is normally degenerate or otherwise pathological,
so ProMesh cannot determine its orientation reliably.

Possible causes include:

    * coincident vertices
    * zero or nearly zero volume
    * collapsed edges or faces
    * invalid volume connectivity
    * severely distorted geometry

For a healthy volume mesh, this number should normally be zero.
]])

ClearSelection(mesh)

local unorientableCount = SelectUnorientableVolumes(mesh)

printCount(
    "Selected unorientable volumes",
    unorientableCount,
    totalVolumes
)

if unorientableCount == 0 then
    print("Result: no unorientable volumes were detected.")
else
    print("Result: these volumes are pathological and require inspection.")
end

if saveSelections then
    local unorientableFile =
        outputPrefix .. "_selected_unorientable.ugx"

    if SaveMesh(mesh, unorientableFile) then
        print("Saved unorientable selection to: " .. unorientableFile)
    else
        print("WARNING: could not save: " .. unorientableFile)
    end
end

ClearSelection(mesh)

-- ============================================================================
-- 5. Orientation test on a cloned mesh
-- ============================================================================

separator("5. INCORRECT VOLUME ORIENTATION")

print([[
This orientation test is performed on a cloned copy of the loaded mesh.
The original input mesh is not modified.

FixVolumeOrientation checks the selected volumes and reverses the vertex
ordering of volumes that have incorrect orientation.

The returned number means:

    0 = all tested volumes were already oriented correctly

   >0 = that many volumes had incorrect orientation and were corrected
        in the cloned mesh

Unorientable volumes and incorrectly oriented volumes are different:

    Unorientable:
        ProMesh cannot determine a reliable orientation.

    Incorrectly oriented:
        ProMesh can determine the orientation, but the vertex order is
        reversed and can be corrected.
]])

-- The original mesh is copied before modifying orientation.
local orientationCopy = CloneMesh(mesh)

ClearSelection(orientationCopy)
SelectAllVolumes(orientationCopy)

local incorrectlyOrientedCount =
    FixVolumeOrientation(orientationCopy)

printCount(
    "Incorrectly oriented volumes",
    incorrectlyOrientedCount,
    totalVolumes
)

if incorrectlyOrientedCount == 0 then
    print("Result: all orientable volumes were already oriented correctly.")
else
    print(
        "Result: incorrectly oriented volumes were corrected " ..
        "in the cloned mesh."
    )
end

if saveSelections then
    local correctedFile =
        outputPrefix .. "_orientation_corrected_copy.ugx"

    if SaveMesh(orientationCopy, correctedFile) then
        print("Saved orientation-corrected copy to: " .. correctedFile)
    else
        print("WARNING: could not save: " .. correctedFile)
    end
end

-- ============================================================================
-- Final summary
-- ============================================================================

separator("QUALITY REPORT SUMMARY")

print("Input grid: " .. inputFile)
print("")

print("Volume composition:")
print("  Tetrahedra:  " .. tetrahedronCount)
print("  Pyramids:    " .. pyramidCount)
print("  Prisms:      " .. prismCount)
print("  Octahedra:   " .. octahedronCount)
print("  Hexahedra:   " .. hexahedronCount)
print("  Other:       " .. otherVolumeCount)
print("  All volumes: " .. totalVolumes)

print("")
print("Detected quality conditions:")

printCount(
    "  Sliver tetrahedra",
    sliverCount,
    tetrahedronCount
)

printCount(
    "  Anisotropic volumes",
    anisotropicCount,
    totalVolumes
)

printCount(
    "  Unorientable volumes",
    unorientableCount,
    totalVolumes
)

printCount(
    "  Incorrectly oriented volumes",
    incorrectlyOrientedCount,
    totalVolumes
)

print("")
print("General interpretation:")
print("  Aspect ratio near 1 is better; near 0 is worse.")
print("  Sliver tetrahedra should ideally be absent or very rare.")
print("  Anisotropic volumes require visual and numerical assessment.")
print("  Unorientable volumes should normally be zero.")
print("  Incorrectly oriented volumes should normally be zero.")

if saveSelections then
    print("")
    print("Diagnostic meshes were saved next to the input grid.")
else
    print("")
    print(
        "Use -saveSelections to save the selected diagnostic meshes."
    )
end

separator("QUALITY REPORT COMPLETED")
