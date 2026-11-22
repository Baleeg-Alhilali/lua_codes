# Project instructions for Codex

This file applies to the `Model` repository. Its purpose is to preserve the
verified lessons from earlier Codex sessions without turning old experiments
into current design rules. Read the relevant source and recent validation
evidence before changing a model. A dated report describes its tested snapshot,
not necessarily the present working tree.

## Orientation and evidence

- The repository root is `Model`. `code/transport`, `code/flow`, and
  `code/coupled` contain the corresponding Lua models; `mesh_generation`
  contains SWC-to-UGX generation and validation. `runs`, `Results`, `output`,
  and `tmp` hold generated data. See `README.md` and `FOLDER_STRUCTURE.txt`.
- Prior session rollouts are in `codex_history/sessions`; recovered reports and
  their provenance are indexed in `codex_history/project_history/INDEX.md`.
  There is no confirmed `archived_sessions` directory or archived-status label.
- Treat experimental filenames, old PDFs, and a successful linear solve as
  evidence of only what was actually tested. Check current code, mesh identity,
  parameters, date, validation log, and physical diagnostics before declaring
  a baseline successful. If two histories conflict, follow the later verified
  result for that exact configuration and explain any remaining uncertainty.
- The working tree may contain substantial user changes and untracked outputs.
  Inspect `git status` first; preserve unrelated work and dated snapshots.
  Use a new output prefix for tests. Never overwrite a historical mesh or
  results set merely to rerun a check.

## Mesh generation: retained practices

- For production multigrid, generate a projector-backed **coarse** UGX and
  use global conforming refinement (`crossOnly=0`). Generator-side
  `_refined_N.ugx` files are flat inspection/validation meshes. The solver
  loads the coarse UGX and builds its hierarchy in memory. Do not substitute
  a flat refined file as the GMG input.
- Keep the Lumen midpoint support ring in production. The center-only option
  created unwanted alternating diagonals. Keep `Apical`, `Basolateral`,
  `MembraneWall`, `Inlet`, `Outlet`, `TripleJunction`, and `OuterWall` distinct;
  terminal MembraneWall annuli are not longitudinal Basolateral faces.
- Preserve the post-refinement cleanup sequence in
  `mesh_generation/refine_neurite_projector.lua`: remove detached vertices,
  fix volume orientation, assert zero unorientable volumes, save, then run
  `validate_ugx_mesh.py`. Generation and standalone refinement validate by
  default; do not disable validation for a production mesh. The October 13
  repair resolved a real projected-refinement inversion, while the validator
  also fixed its own parsing of variable-length constrained-face records.
- Validate newly generated straight and curved controls plus relevant nephron
  traces at coarse, level 1, and level 2 when feasible. Require zero inverted
  or degenerate volumes, sound face topology, subset/interface adjacency, and
  appropriate quality checks. A level-1 pass did not predict level-2 safety
  in earlier projector experiments. Reduced straight/curved fixtures passed
  five refinements; this does **not** prove a full nephron can be materialized
  or solved at level 5. See
  `mesh_generation/MESH_REFINEMENT_VALIDATION_RESULTS_2026-10-13.md`.
- Generated UGX coordinates use millimetres; mesh command length inputs use
  micrometres. Check centerline curvature and nonlocal clearance before
  increasing radius or membrane thickness. Scale the centerline or revise the
  geometry if the tube would self-intersect; do not bypass the preflight to
  make TetGen accept impossible geometry. Terminal extension to OuterWall may
  intentionally bend a nominally straight trace; use local face normals in
  flow code. The validator's `--straight-lumen` check applies only to a mesh
  explicitly required to remain straight.
- Consult `mesh_generation/COMMANDS.txt` and the current Lua defaults for an
  actual run. The October 8 history's anisotropy-1 configuration was a
  verified fixed-anisotropy baseline; the current pipeline exposes adaptive
  anisotropy with straight-path target 4 and bend target 1. Neither result
  licenses arbitrary anisotropy changes. Revalidate the exact configuration.

## Historical fixes and experiments to keep separate

- The September 29 progressive projector and partial radial-correction tests
  are historical. Several corrections passed level 1 but inverted or collapsed
  cells at level 2; even the orientation-safe 25% variant did not improve the
  surface enough to retain. Do not restore a displacement cap or independently
  move connected surface vertices as a routine fix.
- The later October geometry work repaired sharp coarse-centerline turns with
  pre-smoothing and a deterministic transported frame, permitting exact
  projected child placement on its tested configuration. Earlier
  exact-projection failures before that repair do not mean exact projection is
  universally unsafe. Recheck the current generator and source snapshot
  before changing projector behavior.
- `crossOnly=1` produced useful diagnostic meshes after connectivity fixes,
  but remained experimental for the coupled GMG workflow. The apparent
  `*_active_level` surface and a valid local topology test did not establish a
  globally conforming solver hierarchy. Do not select it for production
  without a full geometry and discretization validation.
- The old small-square portal, complete-face portal refinement, reduced
  BoxRefs-3 case, and removal of the Lumen midpoint ring were rejected for
  specific topology, cost, or TetGen failures. The legacy history at
  `codex_history/project_history/legacy_promesh/work/MESH_GEOMETRY_PIPELINE_HISTORY_2026-10-08.md`
  records the experiments; it is not a command to reinstate them.

## Flow, Darcy, and transport

- When starting a simulation, inspect only the first few steps to confirm that
  assembly is finite, the solver is converging, and no timestep retries or
  immediate physical failures occur.  After that initial gate, leave the run
  undisturbed until it completes.  For a long run, check its progress no more
  frequently than about every four minutes, then perform the full algebraic,
  physical, and output checks after completion.  Monitor more closely only
  when a specific warning, instability, retry, or other concrete risk makes it
  necessary.

- For new or substantially edited flow Lua scripts, follow
  `code/flow/FLOW_SCRIPT_STRUCTURE.md`: summary, dependencies, parameters,
  discretization/boundaries, mesh/space, solver, time/initial values, output.
  Keep physical and solver settings in the script's parameters block; expose
  only `-numRefs` and `-output` as runtime inputs for these flow scripts.
  Mesh-generation scripts have their own wider command interface.
- A regular Lumen-only companion mesh was used because NavierStokesFV1 rejected
  constrained objects in the whole mesh. Its stationary Stokes setup passed
  serial and MPI tests at the documented refinements. That result applies to
  the tested lumen solver and mesh, not automatically to a full coupled model.
  Do not apply one multigrid recipe to every problem: the lumen Stokes tests
  used rediscretized coarse operators (`RAP=false`), while a later modified
  Stokes-Darcy interface experiment used RAP/Galerkin (`RAP=true`).
- The September `docs/VELOCITY_FIELD_HANDOFF.md` explicitly reports a large
  flux-balance error in its earlier `velocity_field.lua` test. Do not cite its
  algebraic convergence as physical validation of that velocity field. Use
  later source-specific tests and recomputed inlet, outlet, and interface
  fluxes for any new flow claim.
- In an FV1 inflow test, distinguish the analytic boundary data assembled into
  the continuity equation from its nodal P1 representation written to VTU.
  Directly integrating the coarse exported P1 inlet field can differ from the
  assembled analytic flux and is not, by itself, evidence of mass loss. Check
  the assembled-data inlet flux against the solved outlet/interface flux, and
  separately report convergence of the visualization field under refinement.
- Mesh convergence means that selected physical observables (at minimum flow,
  pressure drop, and a velocity-profile metric) stop changing within an
  explicit tolerance across successive refinements. Exact conservation on each
  level and a small algebraic defect are necessary but do not establish mesh
  convergence. Require at least three levels when feasible; do not call a
  two-level trend converged.
- For the current `Long_tube_V3` Stokes control, this workstation's run scope is
  refinement levels 0--2. Do not start refinement 3 here; hand it to the better
  system with the validated Lumen-only companion and convergence checker. The
  r0--r2 sequence is algebraically and physically valid per level but is not
  mesh converged at the 1% criterion, so keep its ParaView auto-open gated.
- The October 25 Stokes-Darcy technical handoff says the monolithic pressure
  trace/FE Darcy formulation is **not fully validated**. Its `numRefs=2`
  algebraic solve and small true component defects do not resolve possible
  overconstraint from pressure equality plus normal traction, FV/FE row
  consistency, mesh-aligned pressure banding, or Darcy flux reconstruction.
  `DarcyVelocityLinker` is a gradient reconstruction, not proof of the
  assembled conservative boundary flux. Verify the interface laws and a
  control-volume/assembled-flux balance on a small benchmark before treating
  this as a physical production model.
- The November 11 assembled-row diagnostic closes the normal interface-flux
  check for `code/flow/nephron_stokes_darcy_clean.lua`. Immediately before
  pressure-trace elimination, the Stokes and Darcy conservation rows cancel to
  0.00847% at MPI R1 and 0.000213% at MPI R2; serial R1 gives 0.0110%.
  `PressureTraceConstraint` now exposes a read-only, explicitly requested
  capture of those two unreduced row sums, and the clean script gates VTK output
  at 0.1%. Do not use the separately reconstructed `DarcyVelocityLinker` and
  nodal Stokes surface integrals as the mixed FV1/FE conservation test. See
  `docs/ASSEMBLED_INTERFACE_FLUX_VALIDATION_2026-11-11.md`.
- The November 9 straight-tube debugging checkpoint rejected both the old
  pressure-trace/normal-penalty model and the pressure-continuity multiplier:
  the former suppresses or destabilizes normal velocity, while the latter
  double-counts Apical mass flux.  Correcting the stale lumen rim from
  `MembraneWall` to `TripleJunction` makes the Qi-style normal-flux/traction/BJS
  candidate assemble through r2, with component defects below 1e-11 and
  machine-accurate lumen balance.  It is still **not physically valid**:
  Darcy exit mismatches are 97.0%, 79.6%, and 99.3% at r0--r2, and exported
  membrane/inter velocities are strongly non-convergent.  Refined Uzawa/ILU
  GMG produces NaNs and element-patch GMG diverges, so r1/r2 used diagnostic
  SuperLU plus residual correction.  See
  `docs/COUPLED_FLOW_DEBUGGING_2026-11-09.md`; do not promote or auto-open these
  files as verified results.
- The target coupled solver must scale to the eventual full-nephron problem.
  Do not use SuperLU (or another global sparse-direct factorization) as the
  refined/production solver. A direct solve is at most a small-mesh diagnostic;
  in the scalable path it may appear only as the tiny GMG coarse-grid solver.
  Use an outer Krylov method with multigrid and block/saddle-point smoothing,
  and require refinement-independent iteration behavior in addition to true
  residual and physical-conservation checks.
- The later October 27–November 3 coupled-flow session found a concrete
  geometry/input error in
  `code/coupled/coupled_transport_and_flow_only_Na_leak.lua`: its default
  `runs/third_trace_new_build.ugx` has a 1 µm Lumen wall, a 0.5 µm internal
  support ring, and a 1.2 µm Membrane outer wall. The script's current
  `inletRadius=5.0e-4 mm` was set from that support ring, after an earlier
  1 µm assumption was also found wrong for the `Inlet` subset. The user's
  earlier intended Lumen radius for that work was 10 µm. The refined flow-only solve
  passed algebraic residual checks but showed almost zero interior Lumen
  velocity and a porous bypass, so it is **not** a validated physical flow
  result. The subsequent permeability estimate relied on the mistaken radius
  and must not be promoted to a design value. Audit the UGX radial rings,
  subset membership, projected geometry, inlet profile, and output velocities
  together before another coupled run. The newer projector-backed
  `runs/third_trace_generalized3.ugx`, used by the clean flow code and the
  November 12 diffusion-only transport check, has a 25 µm outer tube radius
  and Apical radial fraction 0.6, hence a 15 µm Lumen radius. Its structural
  validation does not settle whether 15 µm replaces the earlier 10 µm model
  target; confirm that physical choice before interpreting transport as
  physiological. No level-0-only solve counts as success for the user's
  refined model.
- The November 12 first transport gate used
  `code/transport/nephron_diffusion_validation.lua` on the generalized3 mesh:
  serial and four-rank MPI R1 closed-Lumen urea diffusion to 0.1 s with a
  nonuniform positive initial profile. Both time steps converged in four
  BiCGStab–GMG iterations in each run. The continuous-P1 mass check changed
  by about 1.49e-7 relative, with serial and MPI integrals agreeing to
  roundoff; independent VTK checks found the same finite, contracting
  concentration bounds. `balancer.firstDistLvl=0` is required for this
  distributed hierarchy: delaying distribution to R1 left the GMG base
  effectively rank-local and failed the iterative solve. UG4 reports 26,025
  concentration DoFs serial and 26,015 MPI; retain that accounting difference
  for investigation, without calling the matching physical integrals a
  failure. This validates short diffusion alone on that mesh, not advection,
  membrane exchange, long-time behavior, or a physiological transport run.
  Keep later transport gates separate and record their own conservation
  evidence.
- The optional whole-Lumen Dirichlet test (`-lumenDirichlet true
  -lumenDirichletValue 1`) on the same R1 mesh kept uniform concentration 1
  and lumen mass exactly constant over two steps to 0.1 s. It confirms the
  boundary value is applied and VTK checkpoints are finite. Because every
  Lumen concentration DoF is fixed, zero iterations and zero mass drift do
  not test diffusion propagation. Use the nonuniform closed-Lumen run above
  as the diffusion-operator evidence.
- The September flow/urea report documents working short transport tests but
  **no proven full 0–60 s refined MPI run**. A tiny GMRES smoke test is not
  full-time convergence. For transport changes, check every adaptive stage,
  concentration bounds, checkpoints, compartment exchange, and mass balance.
  Keep solver conclusions tied to the tested flow and transport blocks.
- Prefer changes in this repository's Lua/Python pipeline and the relevant
  neuro_collection sources. UG4 core or NavierStokes plugin changes have a
  much wider effect and require an explicit need and separate regression
  evidence; old experimental core changes are not a default remedy.

## Reporting a new result

Record the exact source revision or hash, mesh and input identity, full
command/parameters, serial or MPI rank count, validation level, element/DoF
counts, quality and interface checks, solver residuals, physical flux or mass
balance, output path, and known limits in a dated note. Say plainly whether
the result is a geometry pass, an algebraic solve, or a physically validated
model. Update a superseded claim when later evidence overturns it.

For a flow implementation that has passed both its recomputed algebraic-defect
checks and the relevant physical validation checks (including flux or mass
balance), write the final visualization output with a short, descriptive
prefix and automatically open every relevant result dataset in ParaView before
handoff. Open each serial `.vtu` result directly. For parallel output, open the
top-level `.pvtu` dataset rather than its individual partition `.vtu` files.
For coupled Stokes--Darcy tests, "every relevant result dataset" explicitly
means opening the Lumen, Membrane, and Inter compartment datasets together;
never open only the Lumen result after a completed test.
Always launch the compartment `.pvtu` files in ParaView, never their `.pvd`
time-series collection files. When no time step is specified, open the final
available `.pvtu` for Lumen, Membrane, and Inter together.
Use clear short dataset names that identify the model, field or compartment,
and refinement when needed; avoid long experimental-history strings. Do not
present or automatically open a merely algebraically converged result as a
verified physical result. If ParaView cannot be launched or a dataset cannot
be loaded, report the exact blocker and provide the paths that should be
opened.

A serial solve may be retained as a debugging and numerical-reference result,
but it is not a perfect or production-ready implementation. To be described as
perfect, the same model must also solve successfully with MPI on multiple
ranks, pass its independently recomputed component-defect and physical checks,
and write valid top-level parallel `.pvtu` datasets. Compare the MPI observables
against the serial reference; merely launching under `mpirun`, reaching a
Krylov stopping flag, or producing partition files does not satisfy this gate.

## GitHub checkpoints for major work

- When a major model change, mesh-pipeline change, or substantial bug fix reaches
  a reproducible checkpoint, update the affected GitHub repositories during
  that work. Do not leave the only working version on this workstation or
  defer the push indefinitely. This instruction applies separately to `Model`
  and to an external repository such as `neuro_collection` when its source
  actually changed.
- Before committing, inspect each repository's branch, remote, and dirty
  files. Stage only the related source, compact inputs, tests, and dated
  validation notes. Preserve unrelated user changes and keep generated meshes,
  large results, and session rollouts out of ordinary source commits.
- Follow this repository's `development` → `main` workflow in `BRANCHES.txt`:
  commit and push a descriptive checkpoint to `origin/development`; after the
  stated validation passes, merge the reviewed change to `main` and push
  `origin/main`. A large fix whose physical or numerical behavior is still
  unresolved belongs on `development` with its limitation documented, not on
  `main` as a validated result. Follow the affected external repository's own
  branch policy rather than assuming it uses the same branches.
- Record the commit hashes, pushed branches, validation evidence, and any
  remaining limitations in the handoff. If a push or merge cannot be completed,
  report the exact blocker and the local commit state so the checkpoint is not
  mistaken for a GitHub backup.
