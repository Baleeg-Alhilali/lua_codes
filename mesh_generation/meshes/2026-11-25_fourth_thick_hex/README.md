# Full fourth trace with thicker Inter hex cells (2026-11-25)

Same trace, settings and subset names as `mesh_generation/meshes/2026-11-24_fourth_pad100/`
(load the coarse UGX, refine in memory). Only the structured hex shell around
the Membrane changes: 2 layers instead of 4 (the minimum the generator allows).
Built with UG4 `neurites_from_swc.cpp` sha256 300129f9... (linear caps).

| File | Shell / layers | Hex cell radial thickness | Coarse volumes |
|---|---|---|---|
| `fourth_trace_full_pad100.ugx` (old, for reference) | 3 um / 4 | 0.75 um | |
| `fourth_hex2x1p5um.ugx` | 3 um / 2 | 1.5 um | 49,989 |
| `fourth_hex2x2um.ugx` | 4 um / 2 | 2.0 um (thickest valid) | 49,322 |

On the ~41 tight-bend stations the shell is cut to 2.59 um total in both.

Generate (from `mesh_generation`):
`ugshell -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- -swc ../mesh_generation/meshes/2026-11-24_fourth_pad100/fourth_trace_full_pad100_extended_repaired.swc -skipSmoothing 1 -centerlineScale 1 -extendTerminals 0 -radius 10 -membraneThickness 10 -ogrid 12 -BoxRefs 3 -interMesh hex -prismInterThickness <3|4> -interMinLayers 2 -interMaxLayers 2 -numRefs 2 -out <out>.ugx`

Geometry pass only (no solves). Generator validation PASS r0-r2; checker
`check_refinement_quality.lua -mode 2`, 0 inverted/zero volumes r0-r2, min
corner scaled Jacobian r0 / r1 / r2:

| Mesh | tet | hex | prism | pyramid |
|---|---|---|---|---|
| 2 x 1.5 um | 0.0274 / 0.0213 / 0.0211 | 0.914 / 0.954 / 0.952 | 0.49 / 0.469 / 0.470 | 0.245 / 0.081 / 0.073 |
| 2 x 2.0 um | 0.0234 / 0.0038 / 0.0037 | 0.767 / 0.762 / 0.673 | 0.49 / 0.469 / 0.470 | 0.255 / 0.081 / **0.020** |

The 2 um mesh is the thickest valid one but its cap pyramids lose 4x quality
r1 -> r2 (likely inversion near r3) and its worst far-field tet is 6x worse.
The 1.5 um mesh is stable r1 -> r2, like the old 4-layer mesh.
Thicker shells failed: 4.5 um gave 2 inverted cap pyramids at r2; 5, 6, 7 um
segfault in TetGen (5 um reproducibly; also with UG_TETGEN_DOUBLE_PRECISION=1).
Closest coil passes are 56.6 um apart centre to centre (16.6 um between
Membrane surfaces), so shells cannot exceed ~8 um regardless.
Runs: `runs/fourth_thick_hex_2026-11-25/` (not in git).

A prism shell (`-interMesh prism`) was also tried at 3-8 um: all fail at r1/r2
(inverted refined prisms) with far worse coarse tetrahedra; see section 7 of
`docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md`.
