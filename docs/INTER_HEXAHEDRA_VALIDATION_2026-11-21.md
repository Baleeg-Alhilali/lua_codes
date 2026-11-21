# Optional hexahedral near-Inter mesh

Branch: `feature/inter-hexahedra`

The single-nephron mesh pipeline accepts two structured near-Inter topologies:

- `-interMesh prism` retains the previously validated prism implementation and
  remains the default.
- `-interMesh hex` sweeps each radial-angular quadrilateral into one
  hexahedron. The outermost structured layer uses the existing conforming
  pyramid/tetrahedron transition into the far tetrahedral Inter box.

For the hex option, Basolateral and near-Inter refinement is linear. Projecting
only the Basolateral face folded one child of every refined Inter hexahedron;
projecting the complete swept volume folded additional children. Apical
refinement retains the NeuriteProjector. The coarse Basolateral geometry is
unchanged, and all interfaces remain conforming.

## Validated command

```sh
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex /Users/alhilaba/Documents/Promesh_swc/Model/mesh_generation/swc_to_neurite_projector_nephron_final_pipeline.lua \
  -interMesh hex \
  -numRefs 2 \
  -out /Users/alhilaba/Documents/Promesh_swc/Model/runs/third_trace_hex_inter_linear_baso_r2_r10_t10.ugx
```

The other inputs used their 2026-11-21 defaults: lumen radius 10 um, membrane
thickness 10 um, structured near-Inter thickness 20 um, O-grid 10, BoxRefs 4,
straight anisotropy 4, bend anisotropy 1, and padding 200 um.

## Validation result

All coarse, R1, and R2 files passed `validate_ugx_mesh.py`.

| Level | Vertices | Volumes | Minimum hex scaled Jacobian | Result |
|---|---:|---:|---:|---|
| R0 | 45,674 | 137,167 | 0.932378 | PASS |
| R1 | 365,254 | 1,104,436 | 0.195412 | PASS |
| R2 | 2,914,403 | 8,878,088 | 0.195412 | PASS |

R0 Inter volume inventory:

- 14,200 hexahedra in the structured near field
- 3,550 pyramids in the conforming transition
- 105,217 tetrahedra in the transition/far field
- 0 Inter prisms

At R2, the validator reported zero inverted or zero-volume tetrahedra,
hexahedra, and prisms; zero unmatched serialized faces; zero duplicate face
geometries; zero non-manifold faces; and correct Apical and Basolateral
interface adjacency.

This validates mesh geometry and topology through R2 for the tested third
trace. It does not yet validate flow convergence or every admissible SWC.
