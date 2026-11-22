# Two-nephron Stokes versus Navier--Stokes comparison — 2026-11-20

## Configuration

- Navier--Stokes code:
  `code/flow/flow_multiple_nephron_navier_stokes.lua`
- Code SHA-256:
  `49b9b1121454fd5f1844d301c505dfbc1088627dd6f1225045abe8b52fde4f44`
- Stokes reference:
  `code/flow/flow_multiple_nephron_50_50.lua`
- Mesh:
  `runs/third_fourth_nephrons_membranewall_fixed.ugx`
- Refinement: R1
- MPI ranks: 8
- Membrane permeability: `2.35e-12 mm^2`
- Inter permeability: `2.35e-10 mm^2`
- Kinematic viscosity: `1.0 mm^2/s`
- Maximum inlet velocity: `3.0 mm/s`
- Lumen diameter: `0.03 mm`
- Inlet Reynolds number: `0.09`

The nonlinear solve starts from the converged Stokes--Darcy state, activates
the Navier--Stokes convective term, and applies Newton with BiCGStab--RAP-GMG
and Uzawa smoothing.  SuperLU remains restricted to the coarse base level.

## Solver result

- Stokes initialization: 27 BiCGStab iterations, 13.533 s
- Navier--Stokes correction: 3 Newton steps, 49.807 s
- Newton inner BiCGStab iterations: 48, 40, and 21
- Combined solver time: 63.339 s
- Total time: 67.433 s
- Final nonlinear defect: `1.261671e-11`
- Newton absolute stopping tolerance: `1e-10`

The initial `1e-14` nonlinear absolute tolerance was rejected after the MPI
defect reached a repeatable roundoff floor near `1e-11`.  The inner linear
absolute tolerance remains `1e-14`.

## Integrated comparison

| Quantity | Stokes | Navier--Stokes | Change |
|---|---:|---:|---:|
| Nephron 1 outlet fraction | 51.852% | 51.851% | -0.001 percentage point |
| Nephron 1 membrane fraction | 48.148% | 48.149% | +0.001 percentage point |
| Nephron 2 outlet fraction | 48.485% | 48.484% | -0.001 percentage point |
| Nephron 2 membrane fraction | 51.515% | 51.516% | +0.001 percentage point |

Navier--Stokes assembled interface mismatch:

- Nephron 1: `1.17905e-9%`
- Nephron 2: `3.70279e-10%`

## Full-field differences

The PVTU meshes have identical coordinates and ordering.  Differences below
are Navier--Stokes minus Stokes.

| Compartment/field | RMS relative difference | Maximum difference relative to Stokes maximum |
|---|---:|---:|
| Lumen velocity | 0.017194% | 0.148686% |
| Lumen pressure | 0.000842% | 0.001993% |
| Membrane velocity | 0.003342% | 0.030255% |
| Membrane pressure | 0.001869% | 0.021321% |
| Inter velocity | 0.170751% | 3.275140% |
| Inter pressure | 0.006888% | 0.178106% |

The largest local effect is in reconstructed Inter velocity, but its
domain-wide RMS change remains below 0.2%.  At Reynolds number 0.09, inertia
has negligible influence on the 50/50 integrated split.  The Stokes model is
therefore a very accurate and substantially cheaper approximation for this
parameter set.

## Navier--Stokes output

- `Results_final/multi_nephron_navier_stokes_r1_lumen_t0000.pvtu`
- `Results_final/multi_nephron_navier_stokes_r1_membrane_t0000.pvtu`
- `Results_final/multi_nephron_navier_stokes_r1_inter_t0000.pvtu`

This is an R1 comparison, not a mesh-convergence claim.  R2 remains required
before claiming convergence of the Navier--Stokes field differences.
