# Two-nephron 50/50 outlet--membrane calibration — 2026-11-19

## Objective

Increase permeability enough that approximately half of each analytic
parabolic inlet flow leaves through its Outlet and half transfers through its
membrane into the shared Inter domain.

## Model

- Code: `code/flow/flow_multiple_nephron_50_50.lua`
- Mesh: `runs/third_fourth_nephrons_membranewall_fixed.ugx`
- Membrane permeability: `2.35e-12 mm^2` (23.5 times baseline)
- Inter permeability: `2.35e-10 mm^2` (23.5 times baseline)
- BJS permeability: unchanged at `1.0e-13 mm^2`
- Inlet flow reference: `pi * radius^2 * maximumVelocity / 2`
- Solver: MPI 8, R1, BiCGStab + RAP GMG + Uzawa

The BJS parameter was held fixed so the calibration changes normal membrane
exchange and shared-Inter redistribution without simultaneously changing the
tangential slip law.

## Calibration history

- Baseline permeability produced only about 2% membrane transfer.
- 100 times baseline reversed both Outlet flows because the imposed Outlet
  pressure then drove fluid into each lumen while the highly permeable
  membranes drained into Inter.  This was rejected as too aggressive.
- 20 times baseline produced 41.055% and 43.985% membrane transfer.
- 23.5 times baseline gives the final near-50/50 result below.

## Final R1 result

| Nephron | Outlet fraction | Membrane fraction | Interface mismatch |
|---|---:|---:|---:|
| 1 | 51.852% | 48.148% | 4.87612e-6% |
| 2 | 48.485% | 51.515% | 1.21174e-5% |
| Mean | 50.169% | 49.831% | -- |

- BiCGStab iterations: 27
- Solver time: 14.073 s
- Total time: 18.231 s

The two differently shaped nephrons straddle the target while retaining one
common membrane permeability and one common Inter permeability.  Assigning
different permeabilities merely to force exactly 50.000% in each nephron was
not necessary.

## ParaView output

- `Results_final/multi_nephron_50_50_r1_lumen_t0000.pvtu`
- `Results_final/multi_nephron_50_50_r1_membrane_t0000.pvtu`
- `Results_final/multi_nephron_50_50_r1_inter_t0000.pvtu`

## Interpretation limit

The fractions are computed from the analytic inlet integral and the solved
Outlet flux; membrane transfer is the conservative remainder.  The separately
reconstructed nodal Apical flux is not used because it is not the conservative
mixed FV1/FE interface measure.  R2 is still required for a mesh-convergence
claim.
