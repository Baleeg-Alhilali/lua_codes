# Standalone flow and multispecies transport checkpoint (2026-11-22)

## Scope

`code/coupled/nephron_flow_transport_standalone.lua` contains its own
Stokes--Darcy setup, quasi-steady flow coupling, transient Na/K/Cl/urea
advection--diffusion, membrane transport, solver configuration, diagnostics,
and compartment output. It does not load another model script.

The flow inlet is specified by the paper flow rate and converted to a circular
Poiseuille profile:

`u_mean = Q/(pi R^2)` and `u_max = 2 Q/(pi R^2)`.

For `Q = 5.0e-4 mm^3/s` (30 nL/min) and the current mesh inlet radius
`R = 0.010 mm`, `u_mean = 1.591549 mm/s` and `u_max = 3.183099 mm/s`.

## Transport parameters

The paper baseline values are from the `Published models` sheet of
`output/parameters/proximal_tubule_transport_parameters.xlsx`; project-only
diffusion and apical leak assumptions are identified in its
`Current model values` sheet.

| Quantity | Na | K | Cl | Urea |
|---|---:|---:|---:|---:|
| Inlet concentration (mM) | 140 | 4.9 | 113.2 | 5.0 |
| Cell initial concentration (mM) | 19.6 | 138.1 | 16.3 | 4.96 |
| Inter initial concentration (mM) | 140.3 | 4.66 | 112 | 4.91 |
| Fluid diffusion (mm^2/s) | 1.33e-3 | 1.96e-3 | 2.03e-3 | 1.67e-3 |
| Cell diffusion (mm^2/s) | 1.70e-6 | 2.50e-6 | 2.60e-6 | 1.80e-6 |
| Apical passive P (mm/s) | 1.00e-3 | 0 | 8.00e-4 | 3.78e-4 |
| Basolateral passive P (mm/s) | 1.40e-6 | 7.20e-4 | 8.00e-4* | 3.60e-4 |
| Convective reflection coefficient | 0.66 | 0.64 | 0.50 | 0.75 |

`*` Weinstein's cell-lateral passive Cl coefficient is zero. A nonzero
8.00e-4 mm/s project leak assumption is used to satisfy the requested passive
basolateral pathway for every solute. This is a sensitivity parameter, not a
paper-derived basolateral Cl value.

Solvent drag is included on both Apical and Basolateral using

`J_drag = (1-sigma) J_v (c_source+c_target)/2`.

The two-sided FV1 membrane operator uses the signed unit-normal component of
the frozen, converged flow velocity, so one assembled flux is removed from the
source compartment and added to the target compartment. The established UG4
velocity-coupling workaround is retained: flow and transport use separate
approximation spaces, the steady flow is solved first, and the converged
Stokes velocity is passed to transport with `ExplicitGridFunctionVector`.
The Darcy velocities use `ExplicitGridFunctionGradient` of the converged
scaled pressure. This avoids the quadrilateral integration-point evaluation
failure encountered when the flow discretization's native velocity UserData
is attached directly to the transport discretization.

For Basolateral solvent drag, a plugin helper evaluates the converged P1
pressure gradient only on topologically adjacent three-dimensional Membrane
volumes, averages the resulting Darcy velocity at membrane vertices, and
stores it in a membrane-only P1 vector grid function. The boundary operator
then uses ordinary explicit P1 interpolation. This avoids coordinate searches,
prevents selection of an Inter element on the opposite side, and works with
distributed meshes because every evaluation uses local or ghost adjacency.

The available spreadsheet values are whole-epithelium effective convective
reflection coefficients; applying the same value at both explicit membranes
is a provisional assumption until apical- and basolateral-specific values are
identified. `-solventDragScale` defaults to 1 and may be set to 0 for a
drag-disabled control without changing passive leak or pump transport.

The Na/K pump uses `K_Na = 3.516 mM`, `K_K = 0.8584 mM`, and

`J_Na = J_Na,max [Na_i/(Na_i+K_Na)]^3 [K_o/(K_o+K_K)]^2`,

with 3 Na exported for every 2 K imported. The workbook's conversion of
10.8 nmol/(s cm^2) to `1.08e-3 mmol/(mm^2 s)` is dimensionally incorrect:
the correct value is `1.08e-7 mmol/(mm^2 s)`, or `0.108 mM mm/s` for the
concentration units used by the PDE. The transporter therefore receives a
maximum cycle flux of `0.108/3 = 0.036 mM mm/s`.

## Serial smoke test before solvent-drag addition

Command:

```text
mpirun -np 1 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex code/coupled/nephron_flow_transport_standalone.lua \
  -numRefs 0 -endTime 0.01 -dt 0.01 \
  -output tmp/standalone_transport_smoke/model
```

Mesh: `runs/third_trace_third_hex_sim_ready.ugx` (72,769 level-0 elements).
The steady flow BiCGStab--GMG solve converged in two iterations. The exact
assembled Stokes/Darcy interface-row mismatch was `1.3036e-6%`. The exported
coarse P1 surface-flow diagnostic had a 25.4168% lumen imbalance and therefore
does not physically validate the level-0 flow field.

The first backward-Euler transport step converged in two Newton iterations;
its linear solves took two and one BiCGStab--GMG iterations. All Na, K, Cl and
urea point values in the final lumen, membrane and inter VTU files were finite.
Each final `.pvtu` declared all four solutes plus `velocity_mm_s` and referenced
an existing piece file.

This pre-drag test is retained as a baseline.

## Solvent-drag serial smoke test

The drag-enabled level-0 test used the same 0.01-second backward-Euler step
and completed successfully. The flow solve again converged in two iterations.
The exact assembled interface rows were `-1.676619e-5` and
`+1.676619e-5 mm^3/s`; their residual was `-9.51e-15 mm^3/s`. The transport
step converged in two Newton iterations, with one BiCGStab iteration for each
linear update.

The output contains `t0000` and `t0001` `.pvtu` files for Lumen, Membrane and
Inter. Every final file declares Na, K, Cl, Urea and `velocity_mm_s`, references
an existing `.vtu` piece, and all concentration values are finite.

A paired control with `-solventDragScale 0` also completed. At 0.01 s, the
maximum absolute pointwise drag-on versus drag-off concentration differences
(mM) were:

| Compartment | Na | K | Cl | Urea |
|---|---:|---:|---:|---:|
| Lumen | 2.846e-2 | 9.764e-3 | 2.927e-2 | 8.960e-4 |
| Membrane | 10.95 | 19.90 | 13.34 | 0.7235 |
| Inter | 4.396 | 6.687 | 4.499 | 0.2628 |

This confirms that the new term is active and affects every solute. It is a
software smoke test, not a physiological validation: the coarse level-0 flow
still has a 25.4168% reconstructed lumen imbalance, and the exported membrane
velocity magnitude reaches 0.3463 mm/s. A refined parameter-validation run is
required before interpreting the large epithelial concentration response.

The same one-step test also passed with `mpirun -np 8`. Each final Lumen,
Membrane and Inter `.pvtu` references exactly eight rank-local `.vtu` pieces;
ParaView loaded every piece and reported finite Na, K, Cl and Urea arrays.
This specifically validates the MPI-safe basolateral velocity handoff that
replaced the serial-only global element search.

## One-second timestep validation

An eight-rank `endTime=1`, `dt=1` test converged algebraically, but its result
was not admissible. The membrane minima were Na `-122.44 mM`, Cl
`-104.59 mM`, and Urea `-2.160 mM`. Disabling the Na/K pump did not remove the
negative values. Disabling both the pump and solvent drag did remove them
(minimum membrane concentrations remained positive), isolating the problem to
the centered solvent-drag concentration `bar(c)=(c_source+c_target)/2` at this
large timestep. Therefore a 100-second run with `dt=1` was deliberately not
started. The model needs either a smaller internal timestep or a
positivity-preserving/upwind solvent-drag concentration before that run can be
considered valid.
