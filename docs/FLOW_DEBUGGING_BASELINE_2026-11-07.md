# Flow debugging baseline — 2026-11-07

## Scope

This checkpoint separates algebraic convergence from physical conservation for
the existing straight-tube controls. It does not validate the coupled
Stokes–Darcy interface formulation.

UG4 executable:
`/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell` (compiled 2026-11-03).
The accepted refinement-0 through refinement-2 runs were serial. Their short
comparison outputs are `Results_final/stokes_lumen_r0`,
`Results_final/stokes_lumen_r1`, and `Results_final/stokes_lumen_r2`. They pass
the algebraic and per-level physical checks below, but they are not accepted as
mesh-converged results.

## Darcy-only control

Command pattern:

```sh
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex code/flow/flow_straight_darcy_test.lua \
  -numRefs N -output /tmp/promesh_flow_debug_darcy_refN
```

| Refinement | Top-level DoFs | BiCGStab iterations | Recomputed defect |
|---:|---:|---:|---:|
| 0 | 3,812 | 1 | 1.8783e-17 |
| 1 | 27,440 | 27 | 5.9196e-14 |
| 2 | 205,364 | 47 | 1.0094e-13 |

The isolated Darcy FE operator converged through refinement 2. Iteration counts
are not grid independent, but no refinement-born failure was reproduced.

## Stokes-only control

The diagnostic description said Apical was a no-slip wall, but the code applied
both no-normal-stress outflow and a wall-slip/interface condition there. The
control now applies outflow only on Outlet and no-slip wall on Apical. The
retained diagnostic configuration uses the lumen-only Laplace stress form,
four pre/post GMG smoothing sweeps, and ILU damping 0.7. The stopping rule
requires either an absolute defect of 1e-12 or a relative reduction of 1e-12.
A separately recomputed true defect above 1e-11 aborts the run.

Command pattern:

```sh
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex code/flow/debug_stokes_refinement_control.lua \
  -numRefs N -output /tmp/promesh_flow_debug_stokes_final_refN
```

### Flux diagnostic correction

`NavierStokesInflow` contributes the analytic quadratic inlet profile directly
to the FV1 continuity equation, while the VTU contains its nodal P1
interpolation. Integrating the coarse exported P1 field therefore does not
reproduce the assembled boundary flux. This is a normal discretization effect,
not a mass-balance failure.

The correct balance compares the analytic inlet data using FV1 boundary
geometry with the solved outlet flux:

| Refinement | Exact circular flux | Assembled inlet flux | Inlet integration error | Discrete relative balance |
|---:|---:|---:|---:|---:|
| 0 | -4.712389e-4 | -4.949019e-4 | 5.0214% | 1.0597e-12 |
| 1 | -4.712389e-4 | -4.762976e-4 | 1.0735% | 1.4068e-13 |
| 2 | -4.712389e-4 | -4.716507e-4 | 0.08739% | 2.1344e-13 |

All fluxes are in mm^3/s. Both solved levels conserve the assembled flux to
near algebraic tolerance. The changing inlet integration error is expected
until the boundary mesh converges.

### Final solver results through refinement 2

The retained configuration uses the lumen-only Laplace form, ILU damping 0.7,
four pre/post smoothing sweeps, rediscretized coarse operators, and the strict
defect limits stated above.

| Refinement | Lumen cells | Algebraic DoFs | BiCGStab iterations | Recomputed defect | Mean pressure drop (mm^2/s^2) |
|---:|---:|---:|---:|---:|---:|
| 0 | 12,000 | 52,052 | 2 | 2.1209e-14 | 1.379350893157e5 |
| 1 | 96,000 | 296,148 | 33 | 7.8220e-13 | 1.248529001306e5 |
| 2 | 768,000 | 1,936,484 | 50 | 3.8119e-13 | 1.210617686205e5 |

All three runs pass the independently recomputed true-defect assertion and
assembled FV1 physical-flux validation.

### Lumen-only companion validation

`mesh_generation/extract_lumen_flow_mesh.lua` derives
`runs/Long_tube_V3_lumen.ugx` without modifying the 18 MB whole-domain source.
The 2.4 MB companion contains 13,013 vertices and 12,000 Lumen prisms. It
retains Apical, Inlet, Outlet, TripleJunction, and the serialized cylinder
projectors. Standalone validation reports zero inverted/degenerate volumes,
zero unmatched or non-manifold faces, exact subset adjacency, and maximum
Apical radius error `1.7764e-15` micrometres.

The r0, r1, and r2 VTUs produced from this companion are byte-for-byte
identical to the corresponding outputs produced by loading the whole mesh.
Their SHA-256 pairs are respectively:

- r0: `363e0cd57394e8b4fcac364d19c0724b78d1da06cf3e0f4b4bf7de7e6cf7a2a2`
- r1: `a19fb5f958256a293e159220c5cc5453610566d45f2b548a7a5f05a9fb72fdd3`
- r2: `955b142762bd4fd0933e926b36d6845ebbb7730f7f8d4c50ca688088ebb38c3f`

### Mesh-convergence state

`code/flow/check_stokes_mesh_convergence.py` compares nested solutions at
shared physical vertices, removes the optimal constant pressure-gauge offset,
and includes the assembled scalar observables from the Lua logs.

| Pair | Velocity relative L2 | Midplane velocity relative L2 | Centerline axial change | Gauge-free pressure relative L2 | Pressure-drop change | Assembled-flow change |
|---:|---:|---:|---:|---:|---:|---:|
| r0 -> r1 | 26.629% | 26.635% | 26.635% | 10.404% | 10.478% | 3.906% |
| r1 -> r2 | 6.018% | 6.017% | 8.887% | 3.075% | 3.132% | 0.985% |

Only assembled flow meets the initial 1% successive-level tolerance at r2.
Velocity and pressure do not, so mesh convergence is **not established**.

### Refinement-3 handoff

On 2026-11-09, an MPI-4 r3 trial built 6,144,000 Lumen prisms and 13,857,732
algebraic DoFs. Distribution at `firstDistLvl=3` reported 1,536,000 local and
4,608,000 ghost elements per rank. BiCGStab reached iteration 9; the observed
defect history included a reduction from `2.177682e-1` to `1.560568e-3` at
iteration 8 followed by a recoverable-style spike to `6.365191e-2`. The user
requested that r3 be stopped and run on a better system. No r3 VTU/PVTU was
written, and this partial run is not a convergence result.

Handoff SHA-256 values:

```text
efd4bc9752bce1e7d28ba97cdbe9b18846cdd8caca6a60e40e00399e38dbf4fe  runs/Long_tube_V3.ugx
6d8903b964fb89fc1aa52659c961ad341282c0ce38ebba42bcd15a1a1fa36392  runs/Long_tube_V3_lumen.ugx
856f47e0b81f4ea046562f9d7458c98c371e8a0d59ce2192c39b5c56f3319944  code/flow/debug_stokes_refinement_control.lua
832733faf61816a4bc4212394999af0c7bdbdae028d028735515ab3b89e1205e  code/flow/check_stokes_mesh_convergence.py
c9384a5753dbad0a56de4d92f712a52dfb2a8f77f14a253ac1a1cc8b2ff70e3a  mesh_generation/extract_lumen_flow_mesh.lua
```

## Current conclusion

- Darcy refinement is not the immediate convergence blocker in this control.
- The Stokes boundary-condition contradiction was real and is fixed in the
  diagnostic.
- The Stokes linear solve now verifies its true assembled defect.
- The assembled FV1 flux is conservative. Raw VTU P1 boundary-flux integration
  is a different diagnostic and must not be mislabeled as assembled mass loss.
- Refinements 0--2 are verified algebraic/physical diagnostics, but mesh
  convergence is not established. Per the user's workstation scope, do not run
  r3 here. Do not open the r0--r2 files in ParaView as final converged results.
- On 2026-11-09 the user explicitly requested inspection of r2, so
  `Results_final/stokes_lumen_r2_lumen_t0000.vtu` was launched in ParaView
  5.11.2. This inspection exception does not change its non-converged status.

## Resume point

1. Copy the source, extraction script, flow script, convergence checker, and
   Lumen companion to the better system; verify their hashes.
2. Run refinement 3 there, preferably after validating a lower-ghost
   distribution at r2 against the serial byte-identical/reference observables.
3. Require true defect <=1e-11, assembled FV1 flux balance <=1%, and compare r3
   with r2 using the same 1% pressure-drop, assembled-flow, common-vertex
   velocity, midplane-profile, and centerline-velocity criteria.
4. Only after convergence passes, keep short final names and automatically open
   the relevant top-level VTU/PVTU datasets in ParaView. For MPI, open the PVTU,
   never its individual partition VTUs.
