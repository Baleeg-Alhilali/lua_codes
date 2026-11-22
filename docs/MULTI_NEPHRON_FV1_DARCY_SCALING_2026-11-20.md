# Two-nephron FV1 Darcy and pressure-scaling test — 2026-11-20

## Objective

Replace the membrane/interstitium `ConvectionDiffusionFE` Darcy operator by
`ConvectionDiffusionFV1`, while preserving the validated pressure-trace
coupling, and compare algebraic pressure scales `1e4` and `1`.

The mesh element shapes are unchanged. FV1 versus FE changes the discrete
operator and control-volume construction, not whether the UGX contains prisms,
tetrahedra, or other element types.

## Inputs

- FE reference: `code/flow/flow_multiple_nephron_50_50.lua`
- Final FV1 candidate: `code/flow/flow_multiple_nephron_fv1.lua`
- Mesh: `runs/third_fourth_nephrons_membranewall_fixed.ugx`
- Mesh SHA-256:
  `14896c733f7e77e997ecb885bfd93b9d94b86103e9bc5d4df1756fb47c2f29b9`
- FV1 source SHA-256:
  `6a673ddf5bb6fbe9607146aaf94c7c4a4ba15a41a1d926dfcae65d46060244f4`
- MPI ranks: 8
- Refinement: R1
- Solver: GMRES(100), RAP GMG, Uzawa split on `p_stokes`, ILU forward
  smoother, SuperLU on level 0 only

Final command:

```sh
mpirun -np 8 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex code/flow/flow_multiple_nephron_fv1.lua \
  -numRefs 1 -output Results_final/multi_nephron_fv1_r1
```

## Solver development result

The FE-tuned BiCGStab/line-Gauss--Seidel solver did not transfer to FV1:

- FV1, scale `1e4`, line smoother: failed at 400 iterations, final defect
  `8.115044e0`.
- FV1, scale `1`, line smoother: failed at 400 iterations, final defect
  `4.710888e0`.
- FV1, scale `1e4`, ILU forward smoother with BiCGStab: failed at 400
  iterations, final defect `3.917989e-1`.
- Default GMRES used a five-vector restart and stagnated; it was stopped at
  the initial failure gate.
- Explicit GMRES(100) with ILU/GMG converged.

The accepted scale-`1e4` run used two GMRES(100) cycles. Solver time was
`56.857 s`; total time was `61.121 s`. The independently recomputed true
defect was `3.776594086084e-11`. Component defects were:

| Component | Defect |
|---|---:|
| u | 9.107007e-12 |
| v | 7.705738e-12 |
| w | 9.746653e-12 |
| Stokes pressure | 2.138160e-13 |
| Scaled Darcy pressure | 3.448054e-11 |

## Interface conservation and flow split

| Nephron | Outlet fraction | Membrane remainder | Interface mismatch |
|---|---:|---:|---:|
| 1 | 46.363% | 53.637% | 1.52723e-8% |
| 2 | 42.649% | 57.351% | 9.72006e-9% |

These fractions use the analytic inlet reference and solved outlet flux. The
separately reconstructed Apical velocities are not the assembled FV1 flux and
remain unsuitable as the interface-conservation acceptance measure.

## Pressure scale 1 versus 1e4

Changing the pressure scale consistently in the Darcy mobility, OuterWall
value, pressure-trace constraint, traction, and output conversion does not
change the physical problem. The R1 scale-1 and scale-1e4 outputs confirmed
this:

- lumen pressure/velocity RMS relative differences: `2.18e-9` / `7.42e-9`;
- membrane pressure/velocity RMS relative differences: `2.18e-9` / `2.34e-9`;
- inter pressure/velocity RMS relative differences: `1.68e-9` / `8.20e-8`;
- largest reported relative field difference: `1.08e-6` in local Inter
  velocity, at single-precision VTK output scale.

Scale 1 stopped after one GMRES(100) cycle in `32.100 s`, but its absolute
defect was only `1.224076e-6`; it satisfied the relative `1e-8` criterion
because the unscaled Darcy unknown is much larger. Scale `1e4` was retained
because it also passed the explicit absolute true-defect gate of `1e-10`.

## FE versus FV1 at R1

FV1 changes the R1 discrete solution materially:

- nephron-1 outlet fraction: `51.852%` FE versus `46.363%` FV1;
- nephron-2 outlet fraction: `48.485%` FE versus `42.649%` FV1;
- lumen velocity RMS difference: `4.105%`;
- membrane pressure/velocity RMS differences: `14.878%` / `14.543%`;
- Inter pressure/velocity RMS differences: `1.286%` / `16.466%`.

Therefore the FV1 implementation is an algebraically and interface-conservative
R1 candidate, not yet a mesh-converged replacement for FE. R2 is required to
determine whether FE and FV1 approach the same physical solution and whether
the GMRES(100) iteration count scales acceptably.

## Outputs

- `Results_final/multi_nephron_fv1_r1_lumen_t0000.pvtu`
- `Results_final/multi_nephron_fv1_r1_membrane_t0000.pvtu`
- `Results_final/multi_nephron_fv1_r1_inter_t0000.pvtu`

