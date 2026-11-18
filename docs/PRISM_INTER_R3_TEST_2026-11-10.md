# Prismatic Inter r3 flow test — 2026-11-10

## Verdict

The r3 diagnostic **passes its per-level algebraic and lumen-conservation
checks**, and its pressure/velocity oscillation measures improve. It does
**not pass the full production acceptance**: r2-to-r3 pressure drop changes by
7.39%, and BiCGStab/GMG work rises to 747 Krylov iterations. Thus the current
sequence is neither 1%-mesh-converged for all required observables nor backed
by a refinement-independent scalable solver.

## Exact run

Serial, one MPI rank, using the final prism coarse mesh and iterative
BiCGStab/GMG/Uzawa path (no global refined SuperLU):

```text
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -logtofile Results_final/prism_inter_r3.log \
  -ex code/flow/working_flow_coupling_simple_long_tube_pressure_trace_darcy_fe_prism_inter.lua \
  -numRefs 3 -KmemM2 1e-19 -KinterM2 1e-17 -KbjsM2 1e-19 \
  -bjsAlpha 0.2 -outerPressurePa 613.2812 -outletPressurePa 2666.4474 \
  -pressureScale 1e4 -output Results_final/prism_inter_r3
```

The in-memory hierarchy contained 4,232, 33,856, 270,848, and 2,166,784
volume cells at r0--r3. Runtime was 19.3163 minutes.

## Algebra and conservation

- Krylov cycles: 348 + 399 = 747 iterations.
- Final relative true defects: u `4.20e-9`, v `4.13e-9`, w `2.51e-9`,
  Stokes pressure `2.34e-14`, scaled Darcy pressure `8.00e-13`.
- Pressure-trace jump L2: `7.94e-10`.
- Mean Stokes pressure: inlet `2716.090 Pa`, outlet `2666.428 Pa`.
- Outward flows: inlet `-1.058033e-3`, outlet `1.049064e-3`, Apical
  `4.460175e-6 mm3/s`.
- Lumen balance: `-4.508350e-6 mm3/s`, or 0.426% of inlet magnitude.

## r2-to-r3 convergence

| Observable | r2 | r3 | Relative change |
|---|---:|---:|---:|
| Pressure drop | 53.627 Pa | 49.662 Pa | 7.394% |
| Inlet-flow magnitude | 1.051295e-3 | 1.058033e-3 mm3/s | 0.641% |
| Outlet flow | 1.054010e-3 | 1.049064e-3 mm3/s | 0.469% |
| Mean radial membrane velocity | 8.731291e-5 | 8.700947e-5 mm/s | 0.348% |

Flow and radial velocity meet 1%, but pressure drop does not.

## Oscillation diagnostics

- Basolateral circumferential pressure-residual RMS: `0.528 Pa`.
- Circumferential alternating pressure mode RMS: `0.0252 Pa`.
- Middle-60% axial alternating mode: `1.82e-4 Pa`.
- Cell-centered membrane vectors: 235,520, with zero nonfinite components.
- Cell tangential/speed RMS: `5.06e-5` (0.00506%).

These are direct nodal-pressure and `CellData` velocity measurements; no
cell-to-point interpolation was used for the velocity audit.

## Outputs and hashes

- `Results_final/prism_inter_r3_lumen_t0000.vtu`:
  `cf2949259af146338dbecf96cc29864708aa7394f8a91f381bb565744d0a3e34`
- `Results_final/prism_inter_r3_membrane_t0000.vtu`:
  `75958c3147f6fd4d80568619484a74814a90c39bb17dc3602ed2152f89b5e448`
- `Results_final/prism_inter_r3_inter_t0000.vtu`:
  `77c6bdf5d14ca1a3d5b44ac243968ec2cc6096e18b8cff39b9863ca9d5864e7f`
- `Results_final/prism_inter_r3.log`:
  `726407cada0e0e124f3875f497350a72b4e32a2f6d150b28e5be5a2073710062`

ParaView auto-open remains gated because the full acceptance criteria fail,
even though the r3 run itself completed and its per-level conservation passed.
