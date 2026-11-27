# Stokes–Darcy flow on the true-scale full fourth trace — 2026-11-24

Status: **algebraic solve with exact assembled interface conservation and a
converging lumen mass balance; not mesh-converged** (pressure drop still
changes 4.2% from r1 to r2). Not opened in ParaView.

## Inputs

- Mesh: `runs/fourth_trace_full_true_scale/fourth_trace_full.ugx` (coarse,
  projector-backed; refined in memory). Generated from `4th_trace_167.swc` at
  1 um voxels, repaired centreline (min bend 28.3 um, clearance 56.6 um,
  942 um of 1134 um, max move 38.8 um), lumen radius 10 um, membrane 10 um,
  ogrid 12, BoxRefs 3, hex Inter with 3 um shell and outward caps, fade off.
  See `docs/TRUE_SCALE_TRACE_AND_CURVED_BRIDGE_2026-11-23.md`.
- Script: `runs/flow_fourth_full/nephron_flow_fourth_full.lua`, a test copy
  of `code/flow/nephron_stokes_darcy_clean.lua` (SHA-256
  `d9d51900...a839f5`, Model HEAD `bf48df1`) changing only the grid path,
  inlet centre (0.26078106, 0.43055714, 0.15035595) mm and direction (+x),
  output prefix, and adding mean inlet/outlet pressure and DoF statistics.
  Copy SHA-256 `7d1f71a75ad15fe8a66db073f4f9fcb5884a290ad65b4c2d1ce73695781c3838`.
- Parameters unchanged from the clean script: Q = 5e-4 mm^3/s (30 nL/min),
  Poiseuille inlet R = 10 um, nu = 1 mm^2/s, K_mem = 1e-13 mm^2,
  K_inter = 1e-11 mm^2, BJS alpha 0.2, outer 6.132812e5 and outlet
  2.6664474e6 mm^2/s^2, BiCGStab + GMG (RAP, Uzawa/line GS smoother,
  SuperLU on level 0), reduction 1e-8.
- Command (Model root), 8 MPI ranks:
  `mpirun -np 8 ugshell -ex runs/flow_fourth_full/nephron_flow_fourth_full.lua -numRefs N`.
  UG4 binaries as of 2026-11-23 (NeuriteProjector fade extension, opt-in
  TetGen precision; neither affects this solve).

## Results

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Elements | 59,448 | 478,392 | 3,843,984 |
| DoFs (all) | 34,060 | 256,003 | 1,984,321 |
| BiCGStab iterations | 2 (level 0 = direct) | 34 | 84 |
| Average reduction / iteration | — | 0.572 | 0.802 |
| Solver time (s) | 3.5 | 5.9 | 102.6 |
| Total time (s) | 4.0 | 7.0 | 108.9 |
| Assembled Stokes/Darcy interface rows (mm^3/s) | ±2.683566e-6 | ±2.6287e-6 | ±2.6163e-6 |
| Interface row mismatch | 2.5e-7 % | 0.0014 % | 8.5e-5 % |
| Apical reabsorption / Q | 0.537 % | 0.526 % | 0.523 % |
| Outlet flow (nodal, mm^3/s) | 5.111e-4 | 5.015e-4 | 4.979e-4 |
| Closure Q − outlet − apical | −2.76 % | −0.82 % | −0.10 % |
| Nodal inlet integral error vs Q | −11.0 % | −2.85 % | −0.72 % |
| Script "lumen balance" (nodal) | 14.8 % | 3.19 % | 0.40 % |
| Pressure drop (mm^2/s^2) | 1.5037e5 | 1.6862e5 | 1.7598e5 |
| vs Poiseuille 1.8059e5 (L = 1.418 mm) | −16.7 % | −6.6 % | −2.6 % |

- Interface conservation in the assembled rows is exact to round-off at every
  level (gate 0.1%).
- The lumen closure against the prescribed Q converges (−2.8%, −0.8%,
  −0.1%); the residual is the P1 nodal reconstruction of the outlet flux. The
  script's own "lumen balance" uses the nodal inlet integral, whose
  interpolation error falls ~4x per level.
- Pressure drop: changes 10.8% (r0→r1) and 4.2% (r1→r2), observed order 1.3;
  Richardson extrapolation 1.8097e5, 0.2% from the Poiseuille estimate. Not
  converged at the 1% criterion with three levels.
- Solver: iteration count grows with refinement (34 → 84; rate 0.57 → 0.80),
  as previously seen on the 24x mesh; the GMG/Uzawa preconditioner is not
  h-robust.

## Outputs

`runs/flow_fourth_full/fourth_full_flow_r{0,1,2}_{lumen,membrane,inter}_t0000.pvtu`
(8 pieces each), logs `r0.log`, `r1.log`, `r2.log`.

## Limits

- No serial reference run for comparison with MPI.
- Not mesh-converged; r3 would be needed (not run on this workstation).
- Physical model unchanged (hydrostatic reabsorption only, no osmotic term;
  BJS imposed as penalized no-slip as noted on 2026-11-23).

## r3 feasibility on this workstation

r2 rerun with memory sampling (8 ranks): peak total resident memory
**11.6 GB** (at output/diagnostics), solver time 103.6 s, 84 iterations
(identical to the first r2 run). r3 has 8x the elements (~30.8 M) and DoFs
(~15.9 M); memory scales with the finest level, so r3 needs about
**90-100 GB**, against 36 GB on this machine. r3 was therefore **not run**.

Estimate for the larger system (same script, 8 ranks): ~90-100 GB total
memory; if iterations keep growing (34 -> 84 -> ~200), solve time on the order
of 30-40 minutes. Re-check the GMG/Uzawa iteration growth before relying on
this; a refinement-independent preconditioner would cut it substantially.

## Membrane permeability for physiological reabsorption (r0 only)

Target: a rat proximal tubule of ~10 mm reabsorbs ~60% of 30 nL/min. This
lumen (1.418 mm including wall tails, all with Apical membrane) should
reabsorb ~8.5% (uniform per length) to ~12% (load-proportional); target ~10%.
No osmotic/oncotic driving force is modelled, so the permeability found here
is an effective hydrostatic value. Sweep with `-membranePermeability`
(runtime parameter of the test script), all other parameters unchanged,
outputs in `runs/flow_fourth_full/apical_sweep_r0/`:

| K_mem (mm^2) | Reabsorbed (assembled rows, mm^3/s) | % of Q | Pressure drop (mm^2/s^2) |
|---:|---:|---:|---:|
| 1e-13 (default) | 2.684e-6 | 0.54 | 1.504e5 |
| 5e-13 | 1.124e-5 | 2.25 | 1.491e5 |
| 1e-12 | 1.881e-5 | 3.76 | 1.480e5 |
| 2e-12 | 2.872e-5 | 5.74 | 1.465e5 |
| 4e-12 | 3.983e-5 | 7.97 | 1.449e5 |
| **7e-12** | **4.879e-5** | **9.76** | 1.435e5 |
| 1e-11 | 5.433e-5 | 10.87 | 1.427e5 |

All solves: 2 iterations (r0), interface-row mismatch < 3e-7 %.

A series-resistance fit, 1/F = a/K_mem + b, gives a ceiling of ~13.5% of Q
for an infinitely permeable membrane: the interstitium (K_inter = 1e-11 mm^2,
~200 um to the OuterWall at 613 Pa) is the limiting resistance. At
K_mem = 7e-12 mm^2 it carries ~74% of the transepithelial resistance, so the
reabsorption is then controlled mainly by the interstitial model rather than
by the membrane. Equivalent membrane hydraulic conductivity
Lp = K/(mu d) = 7e-18 m^2 / (1e-3 Pa s x 1e-5 m) = 7e-10 m/(Pa s), 70x the
default 1e-11 m/(Pa s). r0 only: the r0 -> r2 change of the default case was
-2.5%, so a similar small shift is expected at r1/r2 (not run).

## K_mem = 7e-12 mm^2, original 200 um padding mesh, r0-r2

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Iterations / solver time | 2 / 3.5 s | 38 / 6.4 s | 97 / 120.7 s |
| Reabsorbed (assembled) | 9.76 % | 9.40 % | 9.28 % |
| Closure Q − outlet − reabsorbed | −2.67 % | −0.17 % | +3.68 % |
| Pressure drop | 1.4355e5 | 1.6058e5 | 1.6483e5 |

r2 rerun at reduction 1e-10 (132 iterations, 204 s) changed the outlet flow
only in the 6th digit: the solve is converged, and the r2 closure gap is a
discretization/flux-reconstruction effect, not solver tolerance. Suspected
cause (unconfirmed): concentrated leakage where the membrane/shell terminal
annuli meet the 613 Pa OuterWall a few um from the lumen.

## Final choice: 100 um padding mesh, K_mem = 7e-12 mm^2

Mesh `runs/fourth_trace_full_short_ext/fourth_trace_full_pad100.ugx` (tails
159 / 122 um, lumen path 1.206 mm; r0-r2 geometry validation PASS). Script
`runs/flow_fourth_full_pad100/nephron_flow_fourth_pad100.lua` (SHA-256
`386ac985...f04824f`): grid path, inlet centre (0.360781065, 0.393057145,
0.16285595) mm, direction +x, K_mem default 7e-12 mm^2. 8 MPI ranks.

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Elements | 52,509 | 422,592 | 3,395,856 |
| DoFs | 30,557 | 229,150 | 1,773,935 |
| BiCGStab iterations (avg reduction) | 2 (direct) | 39 (0.617) | 96 (0.822) |
| Solver time | 3.5 s | 5.8 s | 82.2 s |
| Interface row mismatch | 1.7e-8 % | 2.2e-5 % | 2.0e-5 % |
| Reabsorbed (assembled) | 4.606e-5 (9.21 %) | 4.451e-5 (8.90 %) | 4.401e-5 (8.80 %) |
| Outlet (nodal) | 4.671e-4 | 4.561e-4 | 4.365e-4 |
| Closure Q − outlet − reabsorbed | −2.64 % | −0.13 % | +3.89 % |
| Pressure drop (mm^2/s^2) | 1.2756e5 | 1.3838e5 | 1.4055e5 |

- Reabsorption converges (−3.5%, −1.2%; Richardson 8.75% of Q) and lies in
  the 7-15% physiological target band for this length.
- Pressure drop changes 7.8% then 1.5% (Richardson 1.411e5). Poiseuille for
  1.206 mm without leakage is 1.535e5; with ~8.8% linear leakage ~1.47e5, so
  the extrapolated value is ~4% below that estimate.
- The r2 lumen closure gap (+3.9%) reappears as on the 200 um mesh;
  unresolved (see above). Not mesh-converged at 1% for all observables; not
  opened in ParaView automatically.
- Outputs: `runs/flow_fourth_full_pad100/fourth_pad100_kmem7e-12_r{0,1,2}_{lumen,membrane,inter}_t0000.pvtu`.

## OuterWall pressure sensitivity (100 um padding mesh, K_mem 7e-12 mm^2)

Test option `-outerPressure` added to the pad100 test copy only (default
6.132812e5, unchanged results). r0 sweep (`runs/flow_fourth_full_pad100/outer_sweep_r0/`):

| OuterWall (Pa) | Net reabsorbed (% of Q) |
|---:|---:|
| 613 | +9.21 |
| 1500 | +5.36 |
| 2000 | +3.18 |
| 2400 | +1.45 |
| 2666 (= outlet) | +0.29 |
| 2737 (mid-lumen) | −0.02 |
| 2808 (= inlet) | −0.32 |

Net reabsorption is linear in (mean lumen pressure − OuterWall pressure),
~0.43% of Q per 100 Pa, and vanishes at ~2735 Pa = mean lumen pressure.
With the OuterWall at the mid-lumen pressure the lumen conserves its volume
(outflow upstream, backflow downstream); in this hydrostatic-only model no
membrane permeability then gives net reabsorption — that requires the
osmotic/oncotic terms.

Mid-lumen pressure (2.737e6 mm^2/s^2) at r0-r2 (`outer_mid/`):

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Iterations / solver time | 2 / 3.4 s | 42 / 6.0 s | 99 / 84.3 s |
| Net reabsorbed | −0.018 % | +0.005 % | +0.015 % |
| Closure Q − outlet − reabsorbed | −2.74 % | −0.85 % | −0.30 % |
| Pressure drop | 1.334e5 | 1.452e5 | 1.502e5 (Poiseuille 1.535e5) |

With almost no net transmembrane flux the lumen closure converges (as in the
low-permeability default case), whereas with ~9% net reabsorption it jumped
to +3.9% at r2. The r2 closure gap is therefore tied to large transmembrane
flux, consistent with (not proof of) concentrated leakage near the terminal
OuterWall contact; a per-region flux diagnostic is still needed.

## Diagnosis of the r2 lumen closure gap (100 um mesh, K_mem 7e-12)

Tools: `PressureTraceConstraint::set_interface_defect_dump(prefix)` (new,
read-only, neuro_collection; pre-edit copies and hashes in
`rollback_archive/2026-11-24_interface_flux_dump_pre_edit/`) writes per-node
Stokes/Darcy interface rows; `runs/flow_fourth_full_pad100/flux_diag/section_flux.py`
(pvpython) integrates the Lumen velocity through cross-sections every 25 um.

1. Per node, the Stokes and Darcy interface rows cancel to 2e-11: the
   pressure-trace merge is exact node by node.
2. End leakage is real but small: within 50 um of each end the transmembrane
   flux per length is 2-3x the mid-tube value, an excess of ~3.5e-6 mm^3/s
   (0.7% of Q), nearly unchanged from r1 to r2. It does not explain the gap.
3. The lumen flow decreases along the whole tube faster than the assembled
   transmembrane flux: r1 3.6e-5 vs 3.35e-5 per mm; r2 4.9e-5 vs 3.3e-5 per
   mm. The excess equals the Stokes-side Apical boundary flux
   (u.n of the P1 velocity on Apical faces): 3.2e-6 (r1), 2.05e-5 (r2).
4. Including that term, the lumen closure converges like the low-flux cases:
   Q - outlet - reabsorbed - Apical(u.n) = -2.71 %, -0.78 %, -0.20 % (r0-r2).

Cause: `NavierStokesNoNormalStressOutflow` is applied to "Outlet,Apical"; its
continuity term adds u.n * |face| for Apical faces to the Stokes rows. After
the pressure-trace merge the lumen and membrane control volumes around an
interface node form one union CV whose Apical face is internal, so that term
is an outflow to nowhere: fluid leaves the lumen through it but never enters
the membrane (the Darcy rows do not contain it). It is a mass sink of
4.1% of Q at r2 and grows with refinement as the wall-normal velocity is
resolved; with small transmembrane flux (default K_mem or mid-lumen
OuterWall pressure) it is negligible. The same boundary condition's
tangential viscous flux also double-counts with the BJS friction (noted
2026-11-23).

Candidate fix to test: remove Apical from the NoNormalStressOutflow
condition (keep Outlet), leaving the Apical momentum condition as Darcy
traction (Neumann) + BJS friction + normal impedance, and the mass crossing
implicit in the merged rows. Alternatively add the matching u.n intake to
the Darcy side. Neither is applied yet.

## Fix test: NoNormalStressOutflow on Outlet only

Script `runs/flow_fourth_full_pad100/fix_outlet_only/nephron_flow_fourth_pad100_fix.lua`
= pad100 flux-dump copy with `stokesOpenBoundary:add("Outlet")` instead of
`"Outlet,Apical"`. Apical keeps Darcy traction (Neumann), BJS friction and
normal impedance. 100 um mesh, K_mem 7e-12, OuterWall 613 Pa, 8 ranks.

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Iterations / solver time | 2 / 3.4 s | 40 / 5.8 s | 97 / 79.2 s |
| Interface row mismatch | 8e-9 % | 4.2e-5 % | 8.8e-6 % |
| Reabsorbed (assembled) | 9.21 % | 8.90 % | 8.81 % |
| Closure Q − outlet − reabsorbed (before fix) | −2.64 % | −0.13 % | +3.89 % |
| **Closure (fix)** | **−2.71 %** | **−0.78 %** | **−0.20 %** |
| Pressure drop (fix) | 1.2737e5 | 1.3835e5 | 1.4332e5 |

- Lumen closure now converges as in the low-flux cases; the mass sink is gone.
- r2 cross-section fluxes: lumen flow loss 2.98e-5 per mm (before fix
  4.70e-5) against the assembled mid-tube transmembrane flux 3.32e-5 per mm;
  the remaining ~10% difference in the slope is within the accuracy of the
  P1 section integration.
- Reabsorption unchanged (it was always the water entering the membrane).
- Pressure drop changes 7.9% then 3.5%; Richardson extrapolation 1.474e5,
  matching the leakage-corrected Poiseuille estimate (~1.47e5).
- Solver behaviour unchanged (iterations still grow with refinement).

The fix is applied only to this test copy; `code/flow/nephron_stokes_darcy_clean.lua`
and the coupled/transport scripts that share this boundary setup still use
`"Outlet,Apical"`.

## FV1 Darcy test (100 um mesh, K_mem 7e-12, Outlet-only fix)

Script `runs/flow_fourth_full_pad100/fv1_darcy/nephron_flow_fourth_pad100_fv1.lua`:
fixed copy with membrane and Inter Darcy as `ConvectionDiffusionFV1`, Uzawa
forward smoother ILU (was line Gauss-Seidel) and GMRES(100) outer solver,
following `code/flow/flow_multiple_nephron_fv1.lua` (2026-11-20). 8 ranks.

| | FE r0 | FE r1 | FE r2 | FV1 r0 | FV1 r1 | FV1 r2 |
|---|---:|---:|---:|---:|---:|---:|
| Outer iterations | 2 | 40 | 97 | 2 cycles | 2 cycles | 12 cycles |
| Solver time | 3.4 s | 5.8 s | 79.2 s | 12.2 s | 20.4 s | 944.5 s |
| Reabsorbed | 9.211 % | 8.903 % | 8.810 % | 9.259 % | 8.921 % | 8.817 % |
| Interface mismatch | 8e-9 % | 4.2e-5 % | 8.8e-6 % | 2.9e-9 % | 6.7e-5 % | 8.1e-5 % |
| Lumen closure | −2.71 % | −0.78 % | −0.20 % | −2.71 % | −0.78 % | −0.20 % |
| Pressure drop | 1.2737e5 | 1.3835e5 | 1.4332e5 | 1.2734e5 | 1.3834e5 | 1.4332e5 |

- Physics: FV1 and FE Darcy agree to 0.5% (r0) and 0.08% (r2) in
  reabsorption, identically in lumen closure and pressure drop. FV1 Darcy is a
  reliable discretization here and is locally conservative on the same box
  control volumes as the Stokes FV1 lumen.
- Solver: the ILU-Uzawa GMG with GMRES(100) needed 12 restart cycles at r2
  (about 12x the FE solve time); it degrades far faster with refinement than
  the FE setup and is not usable for r3+ as is. A Darcy-block preconditioner
  suited to FV1 (e.g. block/field-split with an AMG or line smoother for the
  Darcy pressure) would be needed before switching.

## Code-style restyle (2026-11-26)

`code/flow/nephron_flow_fourth_pad100.lua` (FE Darcy, SHA-256 `8fe02886...a9a1ff91`),
`code/flow/nephron_flow_fourth_pad100_fv1darcy.lua` (FV1 Darcy, `b40da63c...a6d526b1`)
and `code/flow/nephron_flow_joint_pad100.lua` (`cfb0d6aa...0bebf783`) were rewritten
in the baleeg-simulation-code-style block layout (initiation, grouped parameters
with units, terminal parameters, functions, discretization with per-equation
boundary conditions, inline `solverDesc`, initial conditions, solve, output).
Physics, solver and diagnostics are unchanged:

- Terminal parameters are now `-grid`, `-numRefs`, `-output`. The test options
  `-membranePermeability`, `-interPermeability`, `-bjsPermeability` and
  `-outerPressure` are gone (plain parameters with the same defaults).
- Default output: `runs/flow_fe_<grid>_<day><Mon>_<hour>h` (FV1 script:
  `flow_fv1_...`).
- The `flowExtraApproximationFunctions` hook was removed (only
  `nephron_stokes_darcy_clean.lua` is loaded by the coupled driver).
- No `assert`; the Model root check prints and exits instead.
- Check: old (f56484d) and restyled scripts run at r0 and r1 on 8 ranks
  (`/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell`) give identical flux,
  pressure, interface-row and true-defect printouts digit for digit and the same
  iteration counts (FE 2/40, FV1 2/2 GMRES cycles, joint 2/36). Test outputs
  deleted.
- The GMG base is still SuperLU on level 0; replacing it belongs to the
  preconditioner work.
