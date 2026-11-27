# Leak transport with block solve per step (2026-11-25/26)

Script: `code/coupled/nephron_flow_transport_standalone.lua`.
Status: **algebraic and discrete-conservation pass at r0 and r1 (8 MPI ranks);
r2 blocked by the flow coarse solve.** Software/transport test stage only: leak
transport across Apical and Basolateral for Na, K, Cl and urea, no pump, no
solvent drag, not physiological.

## Model and solver

- Mesh: `Model/mesh_generation/meshes/2026-11-26_fourth_best_mesh/fourth_trace_best_mesh.ugx`
  (sha256 prefix 4777db46794e095a; 10 um hex Inter shell of two 5 um layers;
  generated data, not in git). Inlet centre (0.338568, 0.400919, 0.163024) mm,
  +x, R = 10 um. New `-grid` option (default this path below `MODEL_ROOT`).
- Flow: Stokes FV1 (Lumen) + Darcy FE (Membrane, Inter), pressure trace,
  K_mem 7e-12 mm^2; BiCGStab + GMG (RAP, Uzawa/line-GS, iterative
  BiCGStab/ILU coarse solve).
- Transport: 12 FV1 fields, FullUpwind with Stokes `velocity_ip()` and Darcy
  `gradient()` velocities (shared approximation space, so the transport uses
  the conservative flow face fluxes), Apical/Basolateral leaks, inlet influx,
  outlet outflow, OuterWall bath Dirichlet.
- Block solve per backward-Euler step (exact for the one-way coupling):
  1. flow block: full Stokes-Darcy system, concentrations held
     (`DirichletBoundary:add(fct, subsets)` keeps current values); solved with
     the flow GMG when its defect exceeds 1e-8 of the steady initial defect;
  2. transport block: flow DoFs held, BiCGStab + GMG (RAP, V(3,3) ILU smoother,
     iterative BiCGStab/ILU coarse solve), reduction 1e-10, defect correction
     until |d| < 1e-10. Jacobian rebuilt only when the flow changes.
- No direct solver anywhere.

## Why not one multigrid for all 17 fields

The monolithic GMG (Uzawa smoother, BiCGStab/Uzawa coarse solve) stalled at r1
on both the old and the best mesh: every coarse solve hit 1000 iterations at a
rate of 0.98-1.00, the flow rows jumped by up to 127x, and the solute rows were
pushed back from ~1e-11 to ~1e-9.

## Root cause found on the way: thin Inter cells

Diffusion+leak tests without flow showed that the transport GMG failure at r1
on `mesh_generation/meshes/2026-11-24_fourth_pad100` came from the Inter hex shell around
Basolateral (0.65 um vs 48 um edges, median edge ratio 0.047), not from the
leak: pure Inter diffusion also diverged under ILU-smoothed GMG (FV1 and FE).
On the best mesh (Inter hex min edge 3.67 um, median ratio 0.29) the same test
needs 1 / 6 / 8-10 iterations at r0 / r1 / r2 with plain ILU.

## Results (8 ranks, dt = 0.1 s)

| | r0 | r1 |
|---|---|---|
| DoFs (17 fields) | 146,890 | 1,078,396 |
| steady flow solve | 2.7 s | 46.3 s |
| interface mismatch / reabsorbed | 0.00036 % / 10.26 % | 0.00015 % / 10.05 % |
| flow block per step | held (defect 8.5e-9) | held (defect 6.5e-7) |
| transport its per step | 1 | 6 (all 5 steps) |
| transport defect per step | ~5e-5 -> ~1e-14 | ~1.8e-5 -> ~2e-15 |
| transport time per step | 0.6 s | 2.3-2.5 s |
| max mass-balance error | 1.7e-8 % | 9.8e-8 % |
| positivity | ok | ok |

r2: the steady flow solve aborts because the GMG coarse solve (BiCGStab/ILU on
the Galerkin coarse Stokes-Darcy operator) does not reach 1e-10 in 2000
iterations. Accepting the inexact coarse iterate made the outer BiCGStab
diverge (46 -> 104 -> 196). The transport block itself is not the limit
(diffusion+leak r2: 8-10 its).

## Limits / next

- r2 needs a working flow coarse solve (e.g. Uzawa-preconditioned coarse
  Krylov, or a deeper coarse hierarchy), without a direct solver.
- Small FV1 undershoots in Inter at coarse levels in diffusion-only tests
  (-6.6e-5 r0, -3.1e-7 r1, none r2); none in the coupled runs.
- Solvent drag and transient flow are not included yet (pump: see below).

## Na/K pump stage (2026-11-26, r0 and r1, 8 ranks)

Transporters: Na apical leak (1e-3 mm/s, project) + basolateral NaKPump only;
K basolateral leak (7.4e-4 mm/s) + pump, no apical leak; Cl and urea leak
across both membranes. NaKPump (neuro_collection): 3 Na out : 2 K in,
Jmax = 0.036 mM mm/s, K_Na = 3.516 mM, K_K = 0.8584 mM (2026-11-22 note).
Nonlinear transport block: chord/Newton defect correction with the reused
Jacobian; rebuild only if a step needs more than 2 corrections (never
triggered up to 1 s).

| | r0 | r1 |
|---|---|---|
| corrections per step (linear its) | 2 (1 + 1) | 2 (6 + 4) |
| transport defect per step | ~4e-5 -> ~2e-8 -> ~5e-11 | ~1.5e-5 -> ~1.2e-8 -> ~3e-11 |
| transport time per step | 0.9-1.0 s | 3.7-3.8 s |
| max mass-balance error | 6.4e-9 % | 3.1e-8 % |
| positivity | ok | ok |

r1 at t = 1 s vs the leak-only run: membrane Na 20.3-21.5 mM (leak only
27.0-27.7), membrane K 129.5-130.2 mM (124.7-125.3), Inter Na up to
141.2 mM next to Basolateral (pumped out), Lumen K max 5.39 mM (7.39; no
apical K leak). Still not physiological: no solvent drag, no apical Na
entry pathways other than the leak, no Cl transporters.

### 30 s run (r1, 8 ranks, dt 0.1, commit a2fa7d2)

Output: `Model/Results/2026-11-26_pump_r1_30s/` (pump_r1_{lumen,membrane,inter},
every 1 s, PVD + PVTU, 1.0 GB, log `run_np8_r1.log`). All 300 steps
converged with 2 corrections (linear its 5-6 + 3-4), no Jacobian rebuild,
mean 3.3 s per step, total 1165 s; max mass-balance error 8.3e-8 %;
positivity ok. Membrane at t = 10 / 20 / 30 s: Na 23.5-26.7 / 25.1-28.5 /
26.1-29.4 mM, K 83.9-85.5 / 64.2-66.0 / 56.9-58.7 mM, Cl 92-94 / 108-111 /
111-115 mM. K keeps draining through the basolateral leak faster than the
pump returns it: without a membrane potential (no electrodiffusion) the cell
cannot hold K, so this is a transport/solver test, not a physiological state.
