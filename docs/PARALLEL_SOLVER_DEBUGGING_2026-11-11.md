# Parallel coupled-flow solver debugging — 2026-11-11

## Status

The projector-backed `third_trace_generalized3.ugx` mesh is retained. The
one-rank Uzawa/line-GMG solve is the numerical reference only. No tested
eight-rank R1 configuration satisfies the algebraic gate, so the coupled-flow
implementation is not yet perfect or production-ready.

The former R1 NaN was traced to a local diagnostic modification in
UG4's distribution layer and removed: `npSurfParams` was explicitly omitted
from global-attachment serialization. Restoring that attachment makes the R1
operator finite on eight ranks. A second, independent bug was then confirmed in
the locally modified Navier--Stokes wall-sliding/BJS Jacobian. Correcting its
missing boundary shape-function factor makes serial and eight-rank R1 converge
in the same 21 total BiCGStab iterations.

`AGENTS.md` now requires a successful multi-rank solve, independently
recomputed component defects, physical checks, serial/MPI agreement, and valid
top-level PVTU files before a flow implementation can be called perfect.

## Reference configuration held fixed

- mesh: `mesh_generation/meshes/legacy_transport_meshes/third_trace_generalized3.ugx`
- refinement: R1 in memory with the serialized projector
- ranks: 8 for MPI probes
- permeability: `Kmem=1e-19 m^2`, `Kinter=1e-17 m^2`, `Kbjs=1e-19 m^2`
- BJS alpha: `0.2`
- pressure scale: `1e4`
- outer pressure: `613.2812 Pa`
- outlet pressure: `2666.4474 Pa`

The serial R1 reference has 698,832 elements and 324,094 algebraic DoFs. It
passes its recomputed component-defect gate, but its physical flow diagnostics
remain invalid (including inward Apical flow and about 5.56% lumen imbalance).

## MPI findings

### Distribution at the finest level (`firstDistLvl=1`)

Assembly stays finite, but every rank retains about 611,458 fine-grid ghosts.
This is not a scalable hierarchy.

- Reference BiCGStab/Uzawa/line GMG: stagnated after repeated 400-iteration
  cycles.
- Explicit agglomerated iterative coarse solver: average reduction
  `1.000229` over the first 200 iterations, then diverged.
- Single-level BiCGStab/Uzawa/line: stable but plateaued at roughly
  `1.4e-2` to `1.8e-2` relative velocity defects after 1,000 iterations.
- Three line sweeps per direction: diverged to velocity defects of order
  `1e6`.
- GMRES(100) with the single-level Uzawa/line preconditioner: preconditioned
  norm remained essentially flat and was stopped.
- Full or block ILU variants generated NaNs; element-patch smoothing diverged.
- Serial refinement followed by final-level partitioning with vertical
  interfaces disabled removed the 611,458 fine ghosts, but produced an invalid
  top-surface view (774,844 surface elements versus 698,832 total hierarchy
  elements, and 363,121 top-surface DoFs versus 324,094 level-1 DoFs). Its
  first BiCGStab step amplified the defect by about `1.97e3`, so this route was
  rejected and is not selectable in the flow script.

### Distribution from level zero (`firstDistLvl=0`)

This produces the desired fine partition (about 87,000 owned fine elements and
zero fine ghosts per rank). Earlier tests produced NaN when the assembled R1
matrix was first applied.

The NaN remains when each of the following is independently removed for a
diagnostic run:

- `PressureTraceConstraint`
- `OneSideP1Constraints`
- BJS terms
- Darcy-pressure traction terms

Those isolation tests correctly excluded the pressure-trace and coupling terms,
but the earlier conclusion that `NavierStokesFV1` itself was defective is
superseded by the attachment audit below.

Additional staged checks establish that the adjusted inlet solution and the
assembled RHS are finite. NaN first appears when the assembled matrix is
applied. A final isolation run removed every boundary, interface, coupling,
and algebraic constraint contribution; the R1 distributed
`NavierStokesFV1` **volume matrix alone** still produced a NaN velocity defect.
This localized the symptom to refined Stokes geometry, but not its cause.

### Confirmed cause of the R1 NaN

The locally modified file
`ugcore/ugbase/lib_grid/parallelization/distribution.cpp` contained:

```cpp
if(name == "npSurfParams") continue; // diagnostic: isolate custom attachment transfer
```

With `firstDistLvl=0`, UG4 distributes the projector-backed coarse grid and
then refines it on the worker ranks. Dropping `npSurfParams` deprives those
ranks of the per-vertex axial, angular, and radial data needed by the
`NeuriteProjector`. The exact pre-edit state is recorded in
`docs/UG4_MPI_PROJECTOR_ATTACHMENT_PRE_EDIT_REFERENCE_2026-11-11.md`.

Removing only this exclusion and rebuilding UG4 changed the same eight-rank R1
test from NaN to finite initial defects: `u=2.532811e2`, `v=2.025405e-2`,
`w=8.686372e2`, `p_stokes=1.283765e-4`, and
`p_darcy_scaled=4.348347e-5`. This confirms a distributed-projector metadata
bug in the local UG4 source state, not a defect in the FV1 flux formula.

### Historical solver failure before the BJS Jacobian fix

Code inspection raised a second hypothesis: `UzawaBase::step` uses block
corrections as inputs to the cross products `B21*u` and `B12*p`. Explicit
unique-to-consistent conversions before those products, followed by conversion
back to unique storage, built successfully but produced an identical R1
trajectory. The experimental edit was therefore rejected and fully reverted;
`uzawa.h` again matches UG4-core `HEAD` with SHA-256
`17a735273e78ab330f17f3b4b5e639fb8123ec01f0b85c8b8a1c99de6423b42d`.
The pre-edit reference is retained in
`docs/UG4_UZAWA_MPI_STORAGE_PRE_EDIT_REFERENCE_2026-11-11.md`.

In the full eight-rank R1 test, the defect was about
`2.24e3` at iteration 148 and still about `1.77e3` at iteration 252; the serial
reference was about `1.93e1` at iteration 148. The run was stopped at iteration
268. Its log is `runs/third_trace_r1_np8_projector_uzawa_fix.log`.

Additional finite-operator probes rejected RAP-off line/Uzawa, point-ILU
Uzawa, element Gauss--Seidel GMG, and overlapping/consistent-interface ILU.
Their 20-iteration average factors were `1.169564`, `1.761237`, `1.078454`, and
`1.368871`, respectively. The experimental overlap-ILU flow profile was
removed.

### Confirmed Navier--Stokes plugin convergence bug

The project-modified
`plugins/NavierStokes/incompressible/fv1/bnd/wall_sliding_fv1.cpp` assembled
the tangential BJS defect from the interpolated boundary velocity,
`sum_sh phi_sh u_sh`, but its Jacobian assigned the full tangential projector
to every nodal column. The Jacobian contribution was missing
`bf->shape(sh)` and therefore was not the derivative of the defect:

```cpp
// incorrect
J(...) += parallelVel(d1,d2) * slidingFactor * boundaryVolume;

// corrected
J(...) += parallelVel(d1,d2) * bf->shape(sh)
          * slidingFactor * boundaryVolume;
```

The exact pre-edit source hash and rollback expression are recorded in
`docs/NAVIER_STOKES_BJS_JACOBIAN_PRE_EDIT_REFERENCE_2026-11-11.md`.

This one-factor correction changed the R1 behavior decisively:

- serial R1: 9 + 12 BiCGStab iterations;
- eight-rank R1: 9 + 12 BiCGStab iterations;
- serial maximum relative component defect: about `2.73e-9`;
- MPI maximum relative component defect: about `2.83e-9`;
- serial/MPI outlet flows: `1.050361e-3` and `1.050358e-3 mm^3/s`;
- serial/MPI lumen balance errors: both about `5.676%`.

The logs are `runs/third_trace_r1_serial_bjs_jacobian_fix.log` and
`runs/third_trace_r1_np8_bjs_jacobian_fix_probe.log`. This confirms the
parallel convergence bug and serial/MPI agreement, but the 5.676% physical
balance error still fails the model-validation gate.

### Controls

- Eight-rank R0 assembled finite and converged in one BiCGStab iteration, with
  component defects below `4e-9` relative. This establishes that MPI execution
  and the coarse coupled operator can work.
- After the projector-attachment fix, serial R0 and eight-rank R0 still converge in one
  iteration. Their outlet flows are `1.416907e-3` and
  `1.416915e-3 mm^3/s`, respectively. Both retain the known physically invalid
  R0 balance and are regression controls only.
- A generic standalone level export was rejected: it changed prism ordering
  and failed the prism Jacobian validation gate. Production GMG must continue
  to load the projector-backed coarse UGX and refine in memory.

## Required next step

The non-finite geometry and R1 parallel-convergence blockers are removed. The
next task is to locate the remaining physical flux-balance discrepancy and
then verify R2 scaling and mesh convergence. Do not promote probe output or
open it in ParaView until both algebraic and physical gates pass.
