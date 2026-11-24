# Two-nephron (third + fourth trace) Stokes–Darcy flow — 2026-11-25

Status: **algebraic solve with exact assembled interface conservation and a
converging lumen mass balance at r0–r2 (8 MPI ranks); not mesh-converged at
the 1% criterion for every observable; no serial reference.** Opened in
ParaView on the user's request for inspection only, not as a validated result.

## Inputs

- Mesh: `inputs/third_fourth_joint_pad100/third_fourth_joint_pad100.ugx`
  (SHA-256 `8afc8636...d29b9e5`, from mesh-generation commit `7407450`,
  copied unchanged into this branch). Coarse projector-backed mesh, refined
  in memory; third trace = nephron 1, fourth trace = nephron 2, per-nephron
  subsets `*_1`, `*_2` and shared `Inter`, `OuterWall`. See that folder's
  README for geometry and validation.
- Script: `code/flow/nephron_flow_joint_pad100.lua` (SHA-256
  `b8391a70...3371a686`), derived from `code/flow/nephron_flow_fourth_pad100.lua`
  on Model `aff59b4`. Changes: subset lists for two nephrons, two Poiseuille
  inlets (R = 10 um; inlet 1 centre (0.173423, 0.432985, 0.237184) mm,
  direction +x; inlet 2 centre (0.528592, 0.637513, 0.188235) mm, direction
  −y; centres are the `Inlet_i` vertex centroids), per-nephron flux/pressure
  diagnostics, true component defects before/after the solve. The runs used a
  test copy differing only in grid/output paths (SHA-256 `2c31c53c...b19d763b5`);
  the committed script reproduces r0 exactly.
- Physics and solver unchanged from the fourth-trace FE script: Q = 5e-4
  mm^3/s (30 nL/min) per nephron, nu = 1 mm^2/s, K_mem = 7e-12 mm^2,
  K_inter = 1e-11 mm^2, BJS alpha 0.2, OuterWall 6.132812e5 and Outlet
  2.6664474e6 mm^2/s^2, NoNormalStressOutflow on Outlets only; FE Darcy;
  BiCGStab + GMG (RAP, Uzawa on p_stokes with line Gauss–Seidel / ILU Schur,
  2+2 smoothing, V-cycle, SuperLU on level 0 only), reduction 1e-8.
- Command (Model root): `mpirun -np 8 <UG4>/bin/ugshell -ex code/flow/nephron_flow_joint_pad100.lua -numRefs N`.

## Results (8 ranks)

| | r0 | r1 | r2 |
|---|---:|---:|---:|
| Elements | 122,420 | 986,224 | 7,930,976 |
| DoFs | 78,196 | 590,295 | 4,586,517 |
| BiCGStab iterations (avg reduction) | 2 (direct) | 36 (0.597) | 101 (0.832) |
| Solver / total time | 19.3 / 21.6 s | 30.4 / 36.0 s | 754 / 786 s |
| Peak resident memory (all ranks) | — | — | 18.1 GB |
| Final true defect u / v / w | ≤1.4e-12 | 4.0e-7 / 5.3e-7 / 7.9e-7 | 1.0e-7 / 8.3e-8 / 5.1e-7 |
| Final true defect p_stokes / p_darcy_scaled | 9.7e-15 / 3.0e-12 | 9.2e-10 / 2.7e-9 | 2.3e-10 / 1.4e-9 |
| Assembled interface mismatch | 5.2e-9 % | 1.1e-3 % | 1.1e-5 % |
| Total reabsorbed (assembled rows) | 8.429e-5 (8.43 % of 2Q) | 8.195e-5 (8.19 %) | 8.125e-5 (8.13 %) |
| Lumen closure 2Q − outlets − reabsorbed | −2.71 % | −0.78 % | −0.20 % |
| Nephron 1 outlet / Q | 94.16 % | 92.45 % | 91.94 % |
| Nephron 2 outlet / Q | 94.41 % | 92.71 % | 92.21 % |
| Nephron 1 pressure drop (mm^2/s^2) | 1.7946e5 | 1.8936e5 | 1.9119e5 |
| Nephron 2 pressure drop (mm^2/s^2) | 1.6497e5 | 1.7715e5 | 1.8121e5 |

- Interface conservation is exact to round-off at every level (gate 0.1 %).
- Lumen closure converges exactly as on the single fourth-trace mesh
  (−2.71 %, −0.78 %, −0.20 %).
- Per-nephron reabsorption is not available from the assembled rows (one
  pressure-trace constraint covers both Apical subsets); Q − outlet gives
  8.06 % (nephron 1) and 7.79 % (nephron 2) of each Q at r2.
- Mesh convergence r1 → r2: total reabsorption 0.85 %, nephron-1 pressure
  drop 0.97 %, nephron-2 pressure drop 2.3 % → **not converged at 1 % for all
  observables**; r3 (~8x memory, ~140 GB) is out of scope on this workstation.
- The fourth trace's pressure drop is higher than on the single-trace mesh
  (1.81e5 vs 1.43e5 at r2) because its joint-repaired lumen path is longer
  (1507 vs 1206 um).
- Solver: iterations grow with refinement (2 → 36 → 101), as for the single
  trace; the GMG/Uzawa preconditioner is still not h-robust. r2 timing was
  taken while other sessions loaded the machine.
- The log's per-nephron "closure Q-outlet-apical(u.n)" line uses the
  Stokes-side P1 u.n on Apical, not the assembled transmembrane flux; use the
  total closure above.

## Outputs

`runs/flow_joint_pad100/` in the main checkout (generated data, not in git):
`joint_r{0,1,2}_{lumen,membrane,inter}_t0000.pvtu` (8 pieces each), per-node
interface dumps `joint_r*_interface_p*.csv`, logs `r0.log`–`r2.log`,
`r2_mem.log`, and `open_joint_r2.py` (ParaView script loading the r2 Lumen,
Membrane and Inter datasets together).

## Preconditioner work in progress (2026-11-24, single fourth trace)

Test harness and logs: `runs/precond_2026-11-24/` (main checkout). Fourth
trace pad100, K_mem 7e-12, Outlet-only fix, 8 ranks, iteration counts only
(machine heavily loaded, times unreliable):

| Darcy / smoother / outer | r1 iterations | r2 iterations | Reabsorbed r2 |
|---|---:|---:|---:|
| FE / Uzawa line GS / BiCGStab (reference) | 40 | 97 | 8.810 % |
| FV1 / Uzawa line GS / BiCGStab | 40 | 92 | 8.818 % |
| FV1 / Uzawa ILU / GMRES(100) (2026-11-24 doc) | 2 cycles (200 inner) | 12 cycles | 8.817 % |
| FV1 / Uzawa ILU / BiCGStab | 141 | — | — |
| FV1 / Uzawa {p_stokes, p_darcy} ILU / BiCGStab | 110 | — | — |
| FE / Uzawa {p_stokes, p_darcy} ILU / BiCGStab | >400 (failed) | — | — |
| FV1 / Uzawa {p_stokes, p_darcy} GS / BiCGStab | 234 | — | — |
| FV1 / plain ILU / BiCGStab | 57 | — | — |

On this mesh FV1 Darcy converges with the FE line-GS solver; the slow FV1
runs came from the ILU/GMRES(100) choice, not from FV1. The r1 → r2 growth
remains; W/F cycles and more smoothing were not yet tested (runs stopped on
the user's request).
