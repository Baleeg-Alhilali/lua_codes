# Two-nephron Stokes--Darcy flow validation — 2026-11-19

## Scope

This checkpoint tests one monolithic Stokes--Darcy solve for the third and
fourth nephron traces inside one shared Inter domain.  The single-nephron
reference script was not modified.

## Repaired mesh

- Production input:
  `runs/third_fourth_nephrons_membranewall_fixed.ugx`
- SHA-256:
  `14896c733f7e77e997ecb885bfd93b9d94b86103e9bc5d4df1756fb47c2f29b9`
- Generator: `mesh_generation/swc_to_neurite_projector_multiple_nephrons.lua`
- Inputs: `3rd_trace_243_third.swc` and `4th_trace_third.swc`
- Geometry: 15 um lumen radius, 10 um membrane thickness, centerline scale 24,
  O-grid 10, graded tetrahedral Inter, projector-backed global refinement.

The repaired coarse mesh contains 40 terminal `MembraneWall_i` faces per
nephron.  The rejected earlier mesh had 3,910 and 4,110 longitudinal faces
incorrectly placed in those subsets.  Coarse validation passed with zero
inverted/degenerate volumes, unmatched faces, duplicate faces, duplicate
subset assignments, and non-manifold faces.

Standalone projected R1 validation also passed.  It contains 160
`MembraneWall_i` faces per nephron, zero unorientable volumes, and the expected
indexed interface adjacency.

## Flow code

- Script: `code/flow/flow_multiple_nephron.lua`
- SHA-256:
  `182de87f19858e605f2dedcede898dc64f829061abb4ea8ebeda0db615a81a01`
- Solver: BiCGStab with RAP GMG, Uzawa smoothing, and SuperLU only on the
  coarse base level.
- MPI ranks: 8
- Each nephron has its own inlet profile, outlet/apical diagnostics, and
  `PressureTraceConstraint`.  The two Darcy membranes share one Inter domain.

## Results

| Level | BiCGStab iterations | Solver time | N1 interface mismatch | N2 interface mismatch |
|---|---:|---:|---:|---:|
| R0 | 2 | 8.677 s | 5.13084e-7% | 3.27116e-7% |
| R1 | 28 | 14.163 s | 8.67745e-6% | 0.000612970% |

R1 flow diagnostics in mm^3/s:

| Nephron | Inlet | Outlet | Reconstructed apical | Nodal lumen balance |
|---|---:|---:|---:|---:|
| 1 | -1.024747e-3 | 1.046641e-3 | 2.870822e-8 | 2.1394% |
| 2 | -1.024747e-3 | 1.044799e-3 | 1.393759e-7 | 1.9704% |

The exact unreduced assembled interface rows cancel independently for both
nephrons and pass the 0.1% output gate.  The separately reconstructed nodal
surface flux is not used as the mixed FV1/FE conservation acceptance test.

## Output

- `Results_final/multi_nephron_flow_r1_lumen_t0000.pvtu`
- `Results_final/multi_nephron_flow_r1_membrane_t0000.pvtu`
- `Results_final/multi_nephron_flow_r1_inter_t0000.pvtu`

All three parallel master datasets were loaded together in ParaView 5.11.2.

## Limitations

This establishes a valid repaired geometry, an MPI algebraic solve, and
per-nephron assembled interface conservation through R1.  It is not yet a
mesh-convergence result: R2 has not been run for this new two-nephron code,
and the reconstructed nodal lumen-balance diagnostic still changes materially
from R0 to R1.
