# Nephron O-grid build snapshot — 2026-10-04

Snapshot name: `O-grid_W_surface_fixed_2026-10-04`

This snapshot preserves the globally conforming multigrid build that removes
the W-shaped Apical/Basolateral surface. It does not replace or modify the
earlier Desktop `O-grid_working_FINAL` snapshot.

Verified production settings:

- anisotropy: 1
- O-grid vertices: 10
- BoxRefs: 4
- uniform/global conforming refinement
- smoothing threshold: 10 degrees
- smoothing strength: 0.85
- smoothing passes: 50
- resampling spacing: 0.5 um
- exact axial and transverse NeuriteProjector placement
- rotation-minimizing frame transport steps: 1024

Validation result:

- coarse: 29,211 vertices and 73,591 volumes
- refinement level 1: 0 indeterminate/degenerate volumes
- refinement level 2: 0 indeterminate/degenerate volumes
- maximum Apical/Basolateral radius residual: 4.15097e-10 mm

The snapshot contains the compiled UG4 source/build tree, active single- and
multiple-nephron Lua pipelines, SWC smoothing and terminal-extension tools,
the detailed investigation notes, validation log, and coarse/L1/L2 inspection
meshes.
