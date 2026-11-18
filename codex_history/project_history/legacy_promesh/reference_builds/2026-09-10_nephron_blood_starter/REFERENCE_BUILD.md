# Nephron + blood starter reference build

Frozen on 2026-09-10 before changing cross-section resolution and cap
classification.

- Outer-box default: `boxRefs = 5`
- Projection selection: `projSH` subset 0 for structured nephron and blood
  tubes; subset 1 for the box and Inter
- Neurite projector: full analytical displacement (no 2% movement cap)
- Verified coarse, refinement 1, and refinement 2 meshes with zero incorrect
  orientations and zero degenerate/unorientable volumes
- Reference mesh: `runs/nephron_blood_full_projector_box5.ugx`

Files in this directory are snapshots and must not be edited while developing
the next version.

Restore into the isolated build with:

    cp neurites_from_swc.cpp /Users/alhilaba/Documents/UG4_promesh/ug4/plugins/neuro_collection/grid_generation/neurites_from_swc.cpp
    cp neurite_projector.cpp /Users/alhilaba/Documents/UG4_promesh/ug4/ugcore/ugbase/lib_grid/refinement/projectors/neurite_projector.cpp
    cmake --build /Users/alhilaba/Documents/UG4_promesh/build -j4

SHA-256 source checksums:

    ea32b7c97a4e089e8eb0e8ba7e580c47594ca392374e65ea5c017ac151ed8172  neurites_from_swc.cpp
    7f7b4799946d2409be3b89f40776f6b823a310a12f17e3ed0fd4f04ffc205132  neurite_projector.cpp
    300dabcf8b7a26aa97fe5d457c614a834addec5ffab75d1b26c2319c2891db44  swc_to_neurite_projector_nephron_blood_pipeline.lua
