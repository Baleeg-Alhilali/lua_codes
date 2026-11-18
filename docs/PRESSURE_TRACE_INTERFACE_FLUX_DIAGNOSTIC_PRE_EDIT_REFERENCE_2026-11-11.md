# Pressure-trace interface-flux diagnostic pre-edit reference — 2026-11-11

This reference records the active `neuro_collection` pressure-trace sources
before adding a read-only diagnostic that captures the separate Stokes and
Darcy interface-row sums immediately before constraint elimination.

The intended change must not alter matrix, right-hand side, defect, transfer,
or solution adjustment. It may only expose values already present in the
unreduced defect when explicitly requested by the flow script.

## Active source hashes

- `pressure_trace_constraint.h`:
  `ed866bcc6e25024575fcad9e699be16e6501749a975b1404f3b9cb9644ccc654`
- `pressure_trace_constraint.cpp`:
  `a781519d08334ccc33904fcbb6d12e183b3592312ea344685404216375aeb7fb`
- `neuro_collection_plugin.cpp`:
  `b91de5744b4306986e8dd2d0586997606d90871fd20449599051acff8f123fda`

## Relevant pre-edit behavior

`PressureTraceConstraint::adjust_defect` first calls `reduce_vector(d, dd)`,
which adds every Darcy trace defect entry to its paired Stokes-pressure entry
and clears the Darcy entry. It then replaces the Darcy entry by the scaled
pressure-continuity defect. No API currently exposes the two contributions
before that reduction, so the assembled normal-flux cancellation cannot be
verified independently from Lua.

