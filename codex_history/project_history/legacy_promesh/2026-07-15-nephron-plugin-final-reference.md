# Nephron mesh plugin final reference

Date: 2026-07-15

## Working locations

- Active plugin: `/Users/alhilaba/Documents/UG4_proMesh/ug4/plugins/neuro_collection`
- Final Desktop backup: `/Users/alhilaba/Desktop/neuro_collection_nephron_final`
- Main implementation: `grid_generation/neurites_from_swc.cpp`
- Header declarations: `grid_generation/neurites_from_swc.h`
- Lua registration: `neuro_collection_plugin.cpp`
- TetGen compatibility change: `test/tetrahedralize_util.cpp`
- Working Lua test workflow: `/Users/alhilaba/Documents/Promesh_swc/swc_to_neurite_projector_membrane_inter_test.lua`
- Automatic smoothing/build pipeline draft: `/Users/alhilaba/Documents/Promesh_swc/swc_to_neurite_projector_nephron_final_pipeline.lua`
- SWC smoother: `/Users/alhilaba/Documents/codes/Thesis/quadratic_interpolate_swc.py`

The Desktop backup was verified to contain the exact working `neurites_from_swc.cpp`. Git metadata and compiled/transient files were excluded. Original version-controlled plugin tests and data were preserved.

## Final domain and boundary layout

- Lumen volume: subset 1, radius 1
- Membrane volume: subset 0, outer radius 1.2 for thickness 0.2
- Inter volume: subset 5
- Basolateral boundary between Membrane and Inter: subset 2
- Apical boundary between Lumen and Membrane: subset 3
- OuterWall box boundary: subset 4
- Box padding: 10 in x, y, and z
- Recommended Inter TetGen quality: 2

The importer exposed to Lua is `import_nephron_with_membrane_and_inter_from_swc`. It includes an `interTetQuality` argument.

## Final mesh architecture

1. Reuse the proven version-2 `import_er_neurites_from_swc` style construction for the conforming Lumen and Membrane hexahedral tube.
2. Keep this version-2 grid untouched so TetGen cannot reorder, remove, or corrupt its vertices.
3. Copy only the nephron exterior surface and padded box surface into an isolated temporary grid.
4. Run TetGen on the temporary grid.
5. Remove temporary tetrahedra that penetrate any Lumen or Membrane cell.
6. Treat Basolateral triangles as an uncrossable topological separator.
7. Flood-fill tetrahedra starting at OuterWall and keep only the exterior component connected to OuterWall without crossing Basolateral.
8. Copy the retained tetrahedra into the original grid as Inter.
9. Use the main/default subset handler for the ProjectionHandler. Assign NeuriteProjector only to Basolateral subset 2 and Apical subset 3.

## Bugs encountered and fixes

### Algebra block-size errors

`AlgebraType("CPU", 1)` and later block size 3 were not compiled in the binary. The working configuration uses an algebra type supported by the installed UG4 build.

### Disconnected SWC / missing soma

The importer reported that SWC points were not connected to a soma. The SWC topology/root data was corrected so the trace formed one valid connected nephron.

### Tube twisting at turns

Fixed and diagnosed through reference-frame warnings and improved SWC smoothing. The strong smoothing workflow uses cubic interpolation and control-point pre-smoothing. The 3rd trace became the clean test trace with no visible flipping or twisting.

### Reference-plane instability

Reference-frame rotation detection was added around the projected reference direction logic. Sudden rotations and sharp generated centerline turns were logged with segment number and coordinates. The projector was restricted to the actual curved tube boundaries rather than volume subsets.

### ProMesh crashes after projector changes

Broad reference-plane/plugin modifications harmed previously correct regions and could crash ProMesh. Those experiments were reverted. The stable approach preserves the original tube generation and makes targeted preprocessing/projector changes.

### Thousands of subsets

An early box/TetGen approach produced thousands of unintended subsets because generated elements inherited or created uncontrolled subset indices. The final implementation uses fixed subset indices 0 through 5 and explicitly names them.

### Inter duplicated inside Lumen and Membrane

The first all-in-one TetGen geometry filled inside the nephron as well as outside. Center-only rejection was insufficient for large tetrahedra. The final approach combines cell penetration sampling with an exact topological flood fill from OuterWall across tetrahedral face adjacency, stopping at Basolateral.

### TetGen corrupted the version-2 tube

TetGen mutates, reorders, inserts, and removes vertices. Reconstructing the tube after TetGen produced distorted geometry. Fixed by running TetGen only in a separate temporary grid and copying the final Inter tetrahedra back.

### Fan-like Inter elements at bends

Low-quality or center-based decompositions sent many edges to one vertex. Using TetGen for the exterior Inter region and `interQuality = 2` produced the accepted mesh quality near turns.

### Membrane fanning after refinement

The ProjectionHandler initially used a separate subset handler that was not reliably preserved after saving/reopening. The final UGX uses the default subset handler and applies the projector only to subsets 2 and 3. Inter, Lumen, Membrane, and OuterWall are not projected as tube surfaces.

### Inter and Membrane appeared to overlap

UGX inspection showed the volume index sets were disjoint; the visual problem was face subset assignment. Faces created implicitly while copying Inter tetrahedra inherited default subset 0, named Membrane. The final fix explicitly assigns all ordinary faces adjacent to Inter tetrahedra to subset 5 while preserving true Basolateral subset 2 and OuterWall subset 4 faces.

The accepted corrected test output was:

- `/Users/alhilaba/Documents/Promesh_swc/runs/third_trace_face_subsets_fixed_ref1.ugx`
- Coarse companion: `/Users/alhilaba/Documents/Promesh_swc/runs/third_trace_face_subsets_fixed_ref1_coarse.ugx`

## Typical validation command

From `/Users/alhilaba/Documents/Promesh_swc`:

```bash
~/Documents/UG4_proMesh/ug4/bin/ugshell \
  -ex swc_to_neurite_projector_membrane_inter_test.lua \
  -swc "./Nephron traces/3rd_trace_243_mesh_ready.swc" \
  -out "./runs/third_trace_face_subsets_fixed_ref1.ugx" \
  -radius 1 \
  -membraneThickness 0.2 \
  -padding 10 \
  -anisotropy 2 \
  -interQuality 2 \
  -numRefs 1
```

## Future multi-SWC work

The current importer supports one SWC filename and one connected nephron geometry. Calling it separately for multiple SWCs creates separate, overlapping boxes rather than one conforming domain. Combining disconnected traces into one SWC can also violate the importer root/soma assumptions.

Proper multi-SWC support should:

1. Accept a list of SWC files.
2. Generate each Lumen/Membrane tube independently.
3. Compute one bounding box around every trace.
4. Insert all Basolateral surfaces as constraints in one temporary TetGen grid.
5. Generate one Inter component outside all membranes and inside OuterWall.
6. Decide whether each nephron uses shared or separately numbered Lumen, Membrane, Apical, and Basolateral subsets.

Do not attempt multi-SWC support by repeatedly calling the existing single-SWC importer into the same box.
