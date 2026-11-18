# MembraneWall subset classification fix — 2026-11-18

## Defect

The common indexed nephron subset pass classified a one-sided structured
Basolateral quadrilateral as `MembraneWall` before checking its preserved
Basolateral identity.  The Inter-side representation uses two constrained
triangles, so ordinary face adjacency sees the structured quadrilateral as
adjacent only to Membrane.  In the former two-nephron output this placed 3,910
and 4,110 faces in `MembraneWall_1` and `MembraneWall_2`, respectively.

Both public generators use this pass: the single-nephron entry point invokes
it with `N=1`, and the multiple-nephron entry point invokes it with the number
of merged traces.

## Repair

`plugins/neuro_collection/grid_generation/neurites_from_swc.cpp` now gives the
preserved Basolateral marker precedence over the general one-sided membrane
rule.  Generation also audits every `MembraneWall[_i]` face before writing the
UGX.  Each such face must belong to the matching Membrane component, must not
touch Lumen, and all of its projector axial parameters must be at one terminal
(`0` or `1`).  Missing terminal wall faces or any longitudinal wall face abort
generation.

## Validation

The rebuilt `neuro_collection` plugin passed these fresh validations:

- Single curved trace, coarse production mesh: 40 `MembraneWall` faces and
  6,990 `Basolateral` faces.  The structural validator reported zero inverted
  or degenerate volumes, zero unmatched or duplicate faces, and
  `VALIDATION: PASS`.
- Two straight traces in one shared Inter box, coarse mesh: 40
  `MembraneWall_i` faces and 960 `Basolateral_i` faces for each nephron.
  `VALIDATION: PASS`.
- The same two-nephron control after one projected refinement: 160
  `MembraneWall_i` faces and 3,840 `Basolateral_i` faces for each nephron,
  consistent with four-way face refinement.  The refined validator reported
  zero inverted or degenerate volumes, zero unmatched or duplicate faces, and
  `VALIDATION: PASS`.

Validation outputs use the prefix
`runs/validation_2026-11-18_*_membranewall_fix` and remain local because
generated meshes are intentionally ignored by Git.

## Separate known limitation

Regenerating the large curved third-plus-fourth-trace case reached the older,
pre-subset Inter construction failure `Could not match a Basolateral quad to
TetGen's two triangles` before subset classification.  That independent
TetGen/interface-matching problem is not caused or hidden by this repair.  The
single curved production case and the focused multi-nephron coarse/refined
regressions above exercise and close the `MembraneWall` classification defect.
