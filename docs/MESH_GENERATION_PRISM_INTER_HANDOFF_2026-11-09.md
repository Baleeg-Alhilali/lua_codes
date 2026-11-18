# Mesh-generation handoff: conforming Basolateral and prismatic Inter

Date: 2026-11-09  
Scope: audit and plan only; no mesh-generator or UG4 source was edited.

## Objective

Investigate the mesh-scale Darcy-pressure oscillation generated inside the
Membrane and most apparent on `Basolateral`. The proposed experiment is to
replace the tetrahedral Inter representation adjacent to the quadrilateral
Membrane surface with a conforming prismatic representation.

Do not change the flow physics while testing topology. Use the same projected
straight tube, physical parameters, solver checks, and r0--r2 comparison.

## Critical finding: source and tested UGX are different snapshots

Tested mesh:

`/Users/alhilaba/Documents/Promesh_swc/Model/runs/straight_nephron_no_bend.ugx`

- Modified: 2026-10-27 12:07:18
- SHA-256: `5f79d67ad7dcfd936c0e7fb0ba12c6a231c7288b889b35e4c73dd5e66ac63894`
- Serialized `NeuriteProjector`, `ProjectionHandler`, and `npSurfParams` exist.

Its coarse volume topology is:

| Region/type | Count |
|---|---:|
| Inter tetrahedra | 14,886 |
| Membrane hexahedra | 460 |
| Lumen hexahedra | 230 |
| Lumen prisms | 230 |
| Pyramids | 0 |
| Total volumes | 15,806 |

This mesh uses the historical nonconforming Basolateral representation: a
Membrane quadrilateral is linked to two Inter triangles by a constraining-face
relationship. The validator accepts `Basolateral` as two one-sided face sets
for this historical topology, and the flow script uses `OneSideP1Constraints`.

The currently compiled generator source is newer:

`/Users/alhilaba/UG4_promesh_ogrid/plugins/neuro_collection/grid_generation/neurites_from_swc.cpp`

- Modified: 2026-11-03 20:35:22
- SHA-256: `098e21d9f47280471e82d94a92b0f38c640d93054f303c419cbf489caf61b034`
- Matching plugin library built: 2026-11-03 20:35:34.

The current C++ common pipeline calls `add_inter_to_version2_nephron` at line
7150. That function describes a different, conforming transition: Membrane
hexes at the Inter handoff are reconstructed with pyramids and tetrahedra so
Basolateral triangles are ordinary shared face objects. The older
`add_isolated_inter_to_version2_nephron` constrained-quad path remains in the
source but is not selected by the current importer.

Therefore the observed oscillation cannot yet be attributed to the current
generator source. First regenerate a new straight mesh with the current build
under a new filename and audit it. Do not overwrite the October 27 mesh.

## Active pipeline

```text
input SWC
  -> quadratic_interpolate_swc.py
  -> unit/radius conversion and geometry preflight
  -> extend_swc_terminals_to_outerwall.py or endpoint-only box fitting
  -> configure_adaptive_nephron_meshing (C++ bridge)
  -> import_nephron_with_membrane_and_inter_from_swc (C++ bridge)
  -> nested Lumen/Membrane O-grid plus NeuriteProjector
  -> OuterWall box and Inter/Basolateral transition
  -> coarse projector-backed UGX
  -> optional flat projected inspection levels
  -> validate_ugx_mesh.py
```

Production GMG loads the coarse projector-backed UGX and refines in memory.
Generator `_refined_N.ugx` files are flat inspection meshes, not GMG inputs.

## Files reviewed

### Model repository

- `mesh_generation/swc_to_neurite_projector_nephron_final_pipeline.lua`
  - Main single-nephron entry point: paths, units, smoothing, geometry
    preflight, OuterWall metadata, adaptive controls, C++ importer, projector
    checks, output validation.
- `mesh_generation/swc_to_neurite_projector_multiple_nephrons.lua`
  - Multi-nephron analogue using the same helper and importer family.
- `mesh_generation/quadratic_interpolate_swc.py`
  - Pre-smooths, interpolates, and resamples centerlines.
- `mesh_generation/extend_swc_terminals_to_outerwall.py`
  - Writes fixed box metadata and extends terminals or fits the box to them.
- `mesh_generation/refine_neurite_projector.lua`
  - Removes detached vertices, fixes orientation, rejects unorientable cells,
    saves, and validates after projected refinement.
- `mesh_generation/promesh_refine_check.lua`
  - Minimal ProMesh GUI-style load/refine/save test.
- `mesh_generation/validate_ugx_mesh.py`
  - Parses tetrahedra, prisms, pyramids, hexahedra, constrained faces,
    topology, subsets, and production interface adjacency.
  - Quality/orientation gates are currently strongest for tetrahedra. A prism
    experiment needs positive-Jacobian/scaled-Jacobian checks for prisms and
    transition pyramids/hexahedra.
- `mesh_generation/Mesh_quality.lua`
  - Element-type, quality, and anisotropy statistics.
- `mesh_generation/COMMANDS.txt` and `mesh_generation/README.txt`
  - Active commands, defaults, units, and output placement.
- `mesh_generation/MESH_REFINEMENT_VALIDATION_HANDOFF.md`
  - Earlier Inter inversion and constrained Basolateral history.
- `mesh_generation/MESH_REFINEMENT_VALIDATION_RESULTS_2026-10-13.md`
  - Orientation-repair fix and r0--r2 structural passes.
- `mesh_generation/tools/remove_detached_outerwall_faces.py`
  - Repair utility, not a normal production stage.
- `mesh_generation/generate_third_trace_ogrid10_from_validated_centerline.lua`
  - Dated wrapper around the same importer.
- `mesh_generation/prepare_long_tube2_flow_control.lua` and
  `mesh_generation/extract_lumen_flow_mesh.lua`
  - Downstream flow-control preparation, not primary generation.

Files under `mesh_generation/archive/` are historical experiments, not current
instructions. Do not restore rejected center-only topology, partial radial
corrections, or cross-sectional-only refinement as production defaults.

### Active UG4/neuro_collection build

- `.../grid_generation/neurites_from_swc.cpp`
  - Adaptive controls: line 123.
  - Structured tube/O-grid creation: line 2019.
  - OuterWall construction: line 3512.
  - Current conforming Inter transition: line 4473.
  - Historical constrained transition: line 4812.
  - Common projector/import/refinement pipeline: line 7015.
  - Public nephron importer: line 7366.
- `.../grid_generation/neurites_from_swc.h`
  - Public declarations and parameter contract.
- `.../neuro_collection_plugin.cpp`
  - Lua bridge registration around lines 1045--1060.
- `.../neuro_collection/CMakeLists.txt`
  - Adds the generator source when not building for VRL.
- `/Users/alhilaba/UG4_promesh_ogrid/tools/tetgen`
  - TetGen dependency for the existing Inter generator.
- UG4 `ugcore/ugbase/lib_grid`
  - Extrusion, tetrahedralization, element creation, projection, refinement,
    orientation, and UGX writing. Avoid core changes unless necessary.

Build facts:

- Build directory: `/Users/alhilaba/UG4_promesh_ogrid/build`
- Runtime: `/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell`
- Release build with ProMesh and neuro_collection enabled.
- The top-level build tree is not one Git repository. The modified
  `neuro_collection` source is not currently versioned as its own repository.
  Create a recoverable patch/snapshot before editing.
- `/Users/alhilaba/Documents/UG4_promesh/ug4/plugins/neuro_collection` is an
  older, non-identical source copy; do not confuse it with the active build.

## Why prisms are reasonable—but not automatically correct

A tetrahedron cannot share a quadrilateral face. A Membrane quad can meet
tetrahedral Inter only through a diagonal split, Membrane reconstruction, or a
nonconforming constraint. The observed pattern follows membrane cells, making
the historical constraint a credible source.

For a straight tube, axially sweeping a triangulated Inter cross-section
creates triangular prisms. Their rectangular lateral faces can match
Basolateral quads one-to-one while preserving axial/circumferential structure.

Acceptance requires:

1. Every Basolateral quad is one ordinary shared face between Membrane and Inter.
2. Both sides share the same four vertex objects—no coincident duplicates.
3. No constrained faces or hanging nodes remain on Basolateral.
4. Prism orientation/Jacobians remain positive through projected r0--r2.
5. Any transition to far-field tetrahedra is conforming.

Possible designs:

- Full swept-prism Inter for the straight control.
- A matched prismatic boundary layer, then conforming pyramids/tetrahedra to
  the OuterWall box. This is more transferable to curved nephrons.

Never connect outer prism quads directly to tetrahedra; use pyramids or a
fully conforming triangulated transition.

## Required work order in the new chat

### Phase 0: reproduce current source before editing

1. Preserve the dirty Model tree and unversioned active C++ source.
2. Record hashes and create a recoverable checkpoint.
3. Rebuild the November 3 source without algorithm changes.
4. Generate `straight_nephron_current_conforming.ugx` under a new name.
5. Generate and validate flat r1/r2 inspection meshes globally.
6. Audit element types, constrained faces, and Basolateral adjacency.
7. Prove whether the current source writes its documented conforming
   pyramid/tetrahedron transition.
8. Run the unchanged flow model through r2 and measure the oscillation.

If the current conforming source removes the oscillation, do not add prisms
yet; retain the simpler validated transition.

### Phase 1: selectable prism prototype

If Phase 0 still oscillates, add an option such as `interMesh = tet | prism`
rather than replacing the established path. Implement the straight case first.

Update:

- `neurites_from_swc.h/.cpp` configuration and construction;
- `neuro_collection_plugin.cpp` Lua bridge;
- production Lua CLI and provenance output;
- prism subsets, projector scope, and refinement inheritance;
- terminal `TripleJunction`, `MembraneWall`, Inlet, Outlet, OuterWall handling;
- validator Jacobian/orientation checks;
- commands and dated validation note.

### Phase 2: controlled A/B validation

At r0, r1, and r2 record:

- element counts by region/type;
- constrained-face counts and Basolateral adjacency;
- minimum quality and Jacobians;
- projected Apical/Basolateral radius error;
- true defects, Krylov iterations, and runtime;
- inlet/outlet/Apical balance;
- pressure drop and outlet-flow convergence;
- Basolateral oscillation amplitude after subtracting the smooth trend;
- identical-range ParaView images.

Measure nodal pressure directly. Inspect cellwise Darcy gradient/velocity
separately and disable ParaView cell-to-point averaging when deciding whether
the pattern is numerical or only visualization.

## Flow baseline to preserve

Script:

`code/flow/working_flow_coupling_simple_long_tube_pressure_trace_darcy_fe_ref2_baseline.lua`

Parameters:

- `KmemM2 = 1e-19`
- `KinterM2 = 1e-17`
- `KbjsM2 = 1e-19`
- `bjsAlpha = 0.2`
- `outerPressurePa = 613.2812`
- `outletPressurePa = 2666.4474`
- `pressureScale = 1e4`

October-mesh projected r2 result:

- mass imbalance: 0.797%;
- pressure-trace jump: `5.028249e-10`;
- true relative defects: approximately `1e-9` to `1e-13`;
- solver cycles: 307 + 149 iterations;
- runtime: 1.72 minutes serially.

The r3 np8 run is not valid for oscillation comparison. Its finest-distributed
mesh assembled (8,092,672 elements; 2,737,041 DoFs), but BiCGStab/GMG diverged
from 64.2 to `2.02e4` and remained near `1e3`. No output was accepted.

## Acceptance gate

- Projector metadata present and used.
- Zero inverted, degenerate, or negative-Jacobian cells through r2.
- No unmatched, duplicate, or non-manifold faces.
- Ordinary two-volume Basolateral adjacency with no constrained pair.
- Correct compartment and terminal subsets.
- r2 mass imbalance below 1% and small verified true defects.
- Oscillation improved or eliminated under a fixed ParaView scale.
- Stable r1-to-r2 physical observables.
- No refined global direct solve; SuperLU only on the tiny GMG coarse level.
- After verification, open Lumen, Membrane, and Inter r2 datasets together
  with short names and summarize solver convergence.

## Suggested opening prompt

Read this handoff completely. Start with Phase 0: reproduce and regenerate the
straight tube using the current November 3 neuro_collection source without
changing its geometry algorithm. Prove which Basolateral transition is written.
Only then decide whether to implement the selectable prismatic Inter prototype.
Preserve dirty trees, historical meshes, projector, subsets, flow parameters,
and validation gates.
