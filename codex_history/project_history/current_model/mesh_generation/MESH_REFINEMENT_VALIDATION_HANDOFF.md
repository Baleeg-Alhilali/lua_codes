# Resolved: straight-tube projected-refinement inversion defect

## Geometry note from flow calibration

The displacement from the original straight SWC centerline is intentional.
The inlet and outlet extensions bend toward centers of subdivided OuterWall
squares to improve terminal mesh quality. Flow/interface code must therefore
use local geometric face normals and must not assume a globally z-aligned
tube. The optional `--straight-lumen` validator is only appropriate for meshes
that are explicitly required to remain cylindrical; it is not a production
validity requirement for this mesh or for nephron meshes.

## Current status (2026-10-13, regenerated mesh)

The mesh-generation defect documented below has been fixed. The regenerated
production files now pass `validate_ugx_mesh.py`:

| File | Volumes | Tetrahedra | Inverted | Zero volume | Result |
|---|---:|---:|---:|---:|---|
| `runs/stright_tube.ugx` | 11,045 | 10,021 | 0 | 0 | PASS |
| `runs/stright_tube_refined_1.ugx` | 88,360 | 80,168 | 0 | 0 | PASS |
| `runs/stright_tube_refined_2.ugx` | 706,880 | 641,344 | 0 | 0 | PASS |

At level 2 there are no unmatched faces, duplicate face geometries,
non-manifold faces, or duplicate subset assignments. All 4,096 `Apical` faces
separate `Lumen` and `Membrane`. The minimum normalized tetrahedron volume is
`0.00127852629478` and the minimum all-pairs edge ratio is `0.070201861`.

The coupled script reconstructed levels 0--2 in memory from the coarse mesh and
completed its serial `numRefs=2` BiCGStab/GMG solve. The hierarchy contained
706,880 top-level volumes and 420,243 top-level DoFs. Final true component
defects were between `1.9e-18` and `4.9e-15`.

The same model also passed with `mpirun -np 2` at `numRefs=1` and `numRefs=2`.
For the two-rank level-2 run, final true component defects were between
`5.95e-18` and `5.77e-15`, and the reported pressure/flow diagnostics matched
the serial run.

Flow-coupling work may continue. The remaining Darcy boundary-velocity output
issue is a flux-reconstruction/postprocessing question, not the old inverted
mesh defect.

## Historical defect (superseded)

### Old status

The flow-coupling work should pause until this mesh-generation defect is fixed.
The source mesh is valid with respect to tetrahedron orientation, but projected
refinement introduces inverted tetrahedra in the `Inter` volume subset.

Mesh under test:

`/Users/alhilaba/Documents/Promesh_swc/Model/runs/stright_tube.ugx`

The refined copies were generated with the same projector-preserving procedure
used by the flow script (`mesh_generation/refine_neurite_projector.lua`).

## Confirmed results

| Level | Total volumes | Tetrahedra | Inverted tetrahedra | Exactly zero-volume tetrahedra |
|---:|---:|---:|---:|---:|
| 0 | 5,620 | 5,108 | 0 | 0 |
| 1 | 44,960 | 40,864 | 1 | 0 |
| 2 | 359,680 | 326,912 | 12 | 0 |

At level 2, all 12 inverted tetrahedra belong to `Inter`. None belong to
`Lumen` or `Membrane`. The smallest absolute signed tetrahedral volume is
`6.94140866863e-10`.

The level-2 inverted tetrahedra are:

| Global volume ID | Signed volume | Center `(x, y, z)` |
|---:|---:|---|
| 24016 | -3.010666430484e-09 | (2.268344747, 3.326098365, 1.408262424) |
| 24020 | -2.609055686966e-09 | (2.267557271, 3.322903336, 1.403701285) |
| 24022 | -2.609055686966e-09 | (2.269420743, 3.321402658, 1.401526819) |
| 274368 | -7.770197463240e-09 | (2.286034585, 3.343338905, 1.491565498) |
| 278272 | -3.855471934228e-09 | (2.261083925, 3.335975378, 1.552116831) |
| 318912 | -1.373056856520e-08 | (2.229208760, 3.296753273, 1.764119704) |
| 318914 | -7.445607737567e-10 | (2.223464365, 3.286580252, 1.770557726) |
| 318915 | -7.445607737577e-10 | (2.231090373, 3.284707878, 1.766352715) |
| 318916 | -6.645273259020e-09 | (2.228320686, 3.292676587, 1.765729174) |
| 318917 | -7.445607737542e-10 | (2.229410920, 3.286715828, 1.766845694) |
| 318918 | -7.445607737549e-10 | (2.226881680, 3.288686383, 1.767896975) |
| 318919 | -6.645273259019e-09 | (2.230849927, 3.290706031, 1.764677893) |

The one level-1 inversion is global volume ID `39864`, signed volume
`-5.956486190057e-09`, centered at
`(2.229015231, 3.289758146, 1.766287448)`. This location develops into the
large level-2 cluster near global IDs 318912--318919.

## Topology checks at level 2

- No unmatched serialized faces.
- No duplicate face geometries.
- No non-manifold volume faces (no face with more than two volume neighbors).
- No duplicate subset assignments.
- `Apical`: 2,048 faces, every face separates `Lumen` and `Membrane`.
- `Inlet`: 128 faces, all adjacent to `Lumen`.
- `Outlet`: 128 faces, all adjacent to `Lumen`.
- `OuterWall`: 15,808 faces, all adjacent to `Inter`.
- `MembraneWall`: 256 faces, all adjacent to `Membrane`.
- `Basolateral` is a constrained/nonconforming interface serialized as two
  one-sided face sets: 4,096 faces on `Inter` and 2,048 on `Membrane`. This is
  not an ordinary shared-face adjacency and must be checked against the
  intended hanging-node construction in the mesh generator.

The direct minimum all-pairs vertex-distance ratio at level 2 is `0.0638074`.
Therefore an earlier ProMesh `SelectAnisotropicVolumes(mesh, 0.05)` report that
selected all volumes must not be interpreted as saying every cell has an edge
ratio below 0.05.

## Why this matters to the solver

Negative-Jacobian tetrahedra can make the local Darcy diffusion contribution
indefinite and can break GMG smoothing/coarse-grid assumptions. The observed
level-2 nonlinear/linear divergence therefore cannot currently be attributed
only to the coupling formulation or solver. Fix the refinement/projection
orientation defect first, regenerate level 2, rerun this validator, and only
then resume the coupled Stokes--Darcy solver tests.

## Reproduce the audit

Reusable validator:

`/Users/alhilaba/Documents/Promesh_swc/Model/mesh_generation/validate_ugx_mesh.py`

Example command:

```bash
cd /Users/alhilaba/Documents/Promesh_swc/Model
python3 mesh_generation/validate_ugx_mesh.py /tmp/stright_tube_ref2.ugx --details
```

The validator exits nonzero when it detects an inverted/degenerate tetrahedron
or a topology/subset-adjacency failure.

## Required acceptance criteria before returning to coupling

1. Zero inverted and zero degenerate volumes at levels 0, 1, and 2.
2. No unmatched, duplicate, or non-manifold faces.
3. Preserve the level-2 `Apical`, inlet, outlet, wall, and volume-subset
   adjacency results listed above.
4. Confirm that the two-sided constrained `Basolateral` representation is
   intentional and that its constraint relations remain valid after refinement.
5. Regenerate the mesh from the corrected generator; do not repair only the
   already-produced refined UGX.
