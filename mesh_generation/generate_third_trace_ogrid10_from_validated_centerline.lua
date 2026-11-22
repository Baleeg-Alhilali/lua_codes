-- Rebuild the validated third-trace nephron with ten circumferential vertices
-- on every O-grid ring.  The centerline and terminal-to-box extension are the
-- same millimetre-scale input used by the validated projector mesh; this file
-- intentionally performs no new smoothing or terminal extension.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

local modelRoot = CurrentWorkingDirectory()
local input = modelRoot .. "/runs/third_trace_new_build_outerwall_extended.swc"
local output = modelRoot
    .. "/runs/ogrid10_2026-11-03/third_trace_ogrid10_projector.ugx"

configure_adaptive_nephron_meshing(
    true, 1.0, 0.02, true, 2.0, 1.5, 6.0, 6.0, false
)

import_nephron_with_membrane_and_inter_from_swc(
    input,
    output,
    1.0 / 1.2, -- lumen radius / outer membrane radius
    4.0,       -- straight-section axial anisotropy
    1,         -- also write a flat level-1 inspection mesh
    0.010,     -- outer-box padding [mm]
    2.0,       -- Inter TetGen quality
    4,         -- outer-box surface refinements
    10,        -- circumferential vertices on every O-grid ring
    false,     -- retain the lumen midpoint support ring
    false      -- refine the complete domain
)

print("O-grid-10 coarse mesh: " .. output)
