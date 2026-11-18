# Assembled Stokes--Darcy interface-flux validation — 2026-11-11

## Outcome

The pressure-trace formulation enforces normal mass conservation in the exact
discrete weak equation. The earlier mismatch came from comparing separately
reconstructed visualization fields, not from the assembled interface rows.

The production diagnostic captures the original Stokes-pressure and stored
Darcy-pressure defect entries immediately before `PressureTraceConstraint`
adds the Darcy row to the Stokes continuity row. It is read-only and is active
only after an explicit `capture_next_interface_defect()` request.

## Validation results

All flux-like defect sums are in the model's consistent mm--s units and have
the scale of `mm^3/s`.

| Run | Stokes row | Darcy row | Sum | Relative mismatch |
|---|---:|---:|---:|---:|
| R1, serial | -1.984714e-05 | 1.984932e-05 | 2.186482e-09 | 0.0110154% |
| R1, MPI 8 | -1.985101e-05 | 1.984933e-05 | -1.680440e-09 | 0.00846526% |
| R2, MPI 8 | -1.957075e-05 | 1.957071e-05 | -4.171442e-11 | 0.000213147% |

The clean flow script rejects a result before VTK output when the relative
assembled interface mismatch exceeds 0.1%.

## Interpretation

The constraint supplies weak normal-flux conservation by retaining the sum of
the two conservation rows while enforcing the scaled pressure trace. It does
not impose pointwise equality between two independently reconstructed velocity
fields. `DarcyVelocityLinker` and the nodal Stokes velocity remain useful for
visualization, but their surface integrals are not the acceptance test for the
mixed FV1/FE interface.

This validation closes the missing algebraic interface-conservation check. It
does not by itself establish mesh convergence of all physical observables or
add the future osmotic-pressure transport law.

## Post-edit source hashes

- `pressure_trace_constraint.h`:
  `444ea53fd8cafc9017ecf7487fe8c9a46c8c5ece925508f1a51767b2402e64ec`
- `pressure_trace_constraint.cpp`:
  `679c3019a990aa8c1210750b5896c625d8605da50fdaa085f223694dac3410a6`
- `neuro_collection_plugin.cpp`:
  `85efaa453f957dbe09bd55d1009de12e4ad857b3d91b7aa5354de2b1af118b62`
- `nephron_stokes_darcy_clean.lua`:
  `89020f16b84023a2d8ebdb0926d8e0744858ebb463644db72343a305c89af50f`

