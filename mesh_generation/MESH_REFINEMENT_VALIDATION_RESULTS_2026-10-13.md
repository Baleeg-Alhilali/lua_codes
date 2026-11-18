# Projected-refinement validation results — 2026-10-13

## Release decision

The level-1/level-2 projected-refinement defect described in
`MESH_REFINEMENT_VALIDATION_HANDOFF.md` is fixed. Newly generated straight,
curved, third-trace, and fourth-trace meshes pass the structural validator at
the coarse level and at global refinement levels 1 and 2.

Mesh generation and standalone refinement now execute validation by default.
A validation failure aborts the Lua script before the mesh can be used for a
solver run.

## Root cause

The standalone ProMesh refinement script called `Refine(mesh)` and immediately
saved the result. Projection changes child coordinates after connectivity is
created. A small number of otherwise non-degenerate Inter tetrahedra therefore
required their vertex ordering to be normalized. The standalone path also left
detached child vertices at constrained terminal interfaces.

The generator's internal refinement path already performed both repairs. The
same sequence is now mandatory in `refine_neurite_projector.lua` after every
level:

1. remove unconnected vertices;
2. select all volumes and run `FixVolumeOrientation`;
3. run `SelectUnorientableVolumes` and assert that the count is zero;
4. save the mesh;
5. run `validate_ugx_mesh.py` and abort on failure.

The obsolete unsafe refinement helper was removed from the generation script.

## Validator correction

The prior Python validator parsed `constrained_triangles` as fixed-width
four-integer records. UGX actually stores the face vertices, a constraint type,
and (unless the type is `-1`) a constraining-object index. That bad stride
shifted later face IDs and produced false coarse-grid topology/interface
errors.

The validator now:

- parses all regular, constraining, and constrained triangle/quad records in
  UGX serialization order;
- understands variable-length constrained-face records;
- includes pyramid and octahedron face layouts;
- checks invalid vertex indices, inverted and exactly-zero tetrahedra,
  serialized-face adjacency, duplicate faces, non-manifold faces, duplicate
  subset assignments, and the production interface contract;
- rejects a scale-free edge ratio below `1e-4` or normalized tetrahedron
  volume `|6V|/Lmax^3` below `1e-10` by default;
- has an explicit `--interface-contract none` mode only for small diagnostic
  fixtures whose terminal subsets intentionally differ from production.

## Regenerated production-style tests

All values below come from newly generated files, not repaired copies of old
refined files.

| Case | Level | Vertices | Total volumes | Inverted / zero tets | Unmatched / duplicate / non-manifold faces | Result |
|---|---:|---:|---:|---:|---:|---|
| straight tube | 0 | 2,842 | 10,335 | 0 / 0 | 0 / 0 / 0 | PASS |
| straight tube | 1 | 19,404 | 82,680 | 0 / 0 | 0 / 0 / 0 | PASS |
| straight tube | 2 | 141,563 | 661,440 | 0 / 0 | 0 / 0 / 0 | PASS |
| curved tube | 0 | 3,303 | 11,848 | 0 / 0 | 0 / 0 / 0 | PASS |
| curved tube | 1 | 23,151 | 94,784 | 0 / 0 | 0 / 0 / 0 | PASS |
| curved tube | 2 | 171,805 | 758,272 | 0 / 0 | 0 / 0 / 0 | PASS |
| third nephron trace | 0 | 10,046 | 41,471 | 0 / 0 | 0 / 0 / 0 | PASS |
| third nephron trace | 1 | 79,508 | 331,768 | 0 / 0 | 0 / 0 / 0 | PASS |
| third nephron trace | 2 | 626,751 | 2,654,144 | 0 / 0 | 0 / 0 / 0 | PASS |
| fourth nephron trace | 0 | 10,725 | 43,264 | 0 / 0 | 0 / 0 / 0 | PASS |
| fourth nephron trace | 1 | 85,107 | 346,112 | 0 / 0 | 0 / 0 / 0 | PASS |
| fourth nephron trace | 2 | 671,741 | 2,768,896 | 0 / 0 | 0 / 0 / 0 | PASS |

Every production-style case also had zero duplicate subset assignments and the
expected Apical, Basolateral, Inlet, Outlet, OuterWall, and MembraneWall
adjacencies.

## Level-5 stress test and scale limit

Reduced straight- and curved-tube projector fixtures were each globally refined
five times. Every refinement step reported zero unorientable/degenerate
volumes. The final straight mesh contained 5,242,880 volumes and the final
curved mesh contained 8,945,664 volumes. Both passed geometry/topology
validation:

- inverted tetrahedra: 0;
- exactly-zero tetrahedra: 0;
- unmatched serialized faces: 0;
- duplicate serialized face geometries: 0;
- non-manifold volume faces: 0;
- duplicate subset assignments: 0.

The level-5 scale-free normalized tetrahedron minima were
`1.28891596597e-06` (straight) and `3.74599556773e-07` (curved), both above the
hard `1e-10` quality gate.

This fixture deliberately keeps terminals inside the box, so it was checked
with the diagnostic interface contract rather than the production OuterWall
terminal contract.

Global refinement multiplies 3-D volume count by approximately eight per
level. Applying five global refinements to the regenerated production third
trace would create about 1.36 billion volumes; the fourth trace would create
about 1.42 billion. Those are not practical inspection files or practical
solver grids on this workstation. Level-5 algorithmic stability must therefore
be regression-tested on reduced straight/curved fixtures, while full nephron
traces are validated at the levels that can actually be materialized. The
per-level in-memory orientation assertion and automatic post-save validator
remain active for any requested production level.

## Saved artifacts

The dated folder `runs/validation_2026-10-13` contains the regenerated coarse,
level-1, and level-2 meshes, generator logs, and a validator report for each
level. Earlier builds were not overwritten.
