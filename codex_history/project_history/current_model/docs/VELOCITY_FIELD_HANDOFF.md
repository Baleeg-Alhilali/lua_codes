# Velocity-field model handoff

Updated on 2026-09-03 after validating the first, lumen-only modelling stage.

## Current requested stage

`velocity_field.lua` now solves only:

- steady incompressible Navier--Stokes in `Lumen`;
- a centerline-aligned parabolic velocity on `Inlet`;
- prescribed outward water leakage on `Apical`;
- the natural zero-traction/no-stress condition on `Outlet`.

There is deliberately no Darcy pressure or Darcy velocity in `Membrane` or
`Inter` at this stage. The earlier coupled implementation is preserved in
`velocity_field_coupled.lua` for a later step.

## Files

- Active stage: `velocity_field.lua`
- Preserved coupled stage: `velocity_field_coupled.lua`
- Original reference, unchanged: `Dif_itterative.lua`
- Mesh: `ProMeshFiles/runs/3rd_trace_third_boxref5.ugx`
- Centerline: `ProMeshFiles/runs/3rd_trace_third_boxref5_outer_radius_1p2.swc`
- Validated parallel output:
  `Results/velocity_field_lumen_after_2_refinements_t0000.pvtu`

All paths above are relative to:

`/Users/alhilaba/Documents/Promesh_swc/Model`

## Physics and boundary conditions

- Stabilized equal-order P1 velocity (`u,v,w`) and P1 pressure (`p_lumen`).
- `AlgebraType("CPU", 4)` groups `u,v,w,p_lumen` at every vertex.
- Constant-density equations divided by density, with kinematic pressure
  `p/rho`, and water kinematic viscosity 6.96e5 micrometre^2/s at 37 C.
- Symmetric-gradient viscous stress (`set_laplace(false)`).
- Small pressure stabilization, default `1e-6`.
- Mean inlet velocity: 800 micrometres/s by default.
- Apical leakage speed: 0.25 micrometres/s by default, directed radially out of
  `Lumen` and ramped smoothly to zero at the inlet and outlet rims.
- `Outlet` is left unconstrained, producing the natural FE zero-traction
  boundary condition.
- The creeping-flow Navier--Stokes limit is used because the Reynolds number is
  approximately 0.0023; the inertial term is negligible at this stage.

The inlet is normalized on the real mesh. A level-dependent compensation for
UG4's P1 inflow trace is applied automatically, giving the requested mean flow
at the default two refinements.

## Solver and refinement

- Default `numRefs = 2`.
- The script registers `projSH` before `LoadDomain`, restoring the two
  UGX-embedded `NeuriteProjector`s on `Basolateral` and `Apical` before any
  refinement is performed.
- The local UG4 loader was corrected in
  `ugcore/ugbase/lib_grid/file_io/file_io_ugx_impl.hpp` so the loaded projection
  handler count is propagated to `LoadDomain`; UG4 was rebuilt afterward.
- Outer solver: distributed BiCGStab.
- Preconditioner: geometric multigrid, F-cycle, RAP matrices.
- Fine-level smoother: overlapping parallel ILU.
- Base level: 0, with only 8,128 scalar DoFs in CPU4 blocks.
- SuperLU is used only on this small level-0 coarse system. The refined
  202,900-DoF system remains distributed and is never gathered for a direct
  factorization.

The earlier P2/P1 parallel failure occurred because full SuperLU gathered all
1,036,576 fine-grid unknowns onto rank 0. The new CPU4/P1-P1 configuration
keeps the refined system distributed and converges without the memory error.

## Projector-refinement validation

Command:

```bash
cd "/Users/alhilaba/Documents/Promesh_swc/Model"
mpirun -np 8 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex code/flow/velocity_field.lua
```

Validated mesh/system size:

- 8 MPI processes;
- 2 refinements;
- 5,100,800 total mesh elements;
- 202,900 lumen Navier--Stokes scalar DoFs in CPU4 blocks;
- 32,448 exported level-2 `Lumen` cells.

The embedded `ProjectionHandler` was verified at runtime. The exported level-2
lumen contains 32,448 hexahedra. A Gauss-point Jacobian check found a minimum
determinant of `9.4253e-4` and no nonpositive determinants, so the projected
lumen mesh contains no detected inverted cells.

BiCGStab converged in 45 iterations to a relative reduction of `1e-12`.
The complete run took approximately 0.50 minutes on the local workstation.

Flux convention: positive is outward from `Lumen`.

```text
Lumen inlet flux:       -1.599600912e+03 micrometre^3/s
Lumen apical leakage:    3.986427528e+02 micrometre^3/s
Lumen outlet flux:      -5.831764181e+01 micrometre^3/s
Lumen balance residual: -1.259275801e+03 micrometre^3/s
```

The inlet mean is `1599.600912 / 1.999606129 = 799.958` micrometres/s. The
projected mesh is valid, but this flow result is not yet quantitatively valid:
the large flux-balance residual persists even with a tighter algebraic solve.
Do not use the projected velocity field for mass-balance conclusions until the
equal-order velocity/pressure discretization is revised and revalidated.

## Literature defaults

1. Human S1 proximal-tubule mean velocity: 48.15 mm/min = 802.5
   micrometres/s, rounded to 800 micrometres/s.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC10117878/
2. Rabbit proximal-tubule water reabsorption: 1.9 nL/(mm min), corresponding
   to approximately 0.24 micrometres/s through a 41.5-micrometre-diameter
   wall, rounded to 0.25 micrometres/s.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC436574/

## Next modelling step

Revise and revalidate the pressure/velocity discretization on the projected
mesh, then review the lumen velocity and leakage output in ParaView. Only after
the flux balance is restored should membrane/interstitial pressure and Darcy
velocity be added back.
