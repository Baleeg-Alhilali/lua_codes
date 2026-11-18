# Straight-tube Stokes--Darcy debugging checkpoint (2026-11-09)

## Scope and status

This checkpoint tests the projector-backed `runs/stright_tube.ugx` through
refinement levels 0--2.  It is an algebraically converged debugging result,
not a physically validated coupled-flow baseline.  Refinement 3 was not run.

The active candidate is
`code/flow/qi_2025_paper_interface_stokes_darcy_simple_long_tube.lua`.  It uses
normal mass conservation, Darcy-pressure normal traction, and tangential BJS
slip.  Trial output was written under `/tmp`; nothing was promoted to
`Results_final`.

## Bugs isolated and changes made

1. The old pressure-trace model sealed `OuterWall`, which forces the steady
   porous velocity toward zero.  Opening it at 0 Pa produced nonzero flow but
   exposed a large conservation error.
2. Its `PressureTraceConstraint` plus a normal Robin penalty either suppressed
   interface velocity or exposed an unstable normal mode.  It was rejected.
3. The pressure-continuity multiplier added its multiplier flux on top of the
   Stokes Apical velocity flux.  The resulting equation double-counted
   interface mass and was rejected.
4. Several experimental scripts incorrectly used `MembraneWall` as the lumen
   rim.  On this mesh the actual rim is `TripleJunction`; the stale map made
   Navier--Stokes boundary assembly touch porous prisms and caused crashes.
5. The Qi-style Darcy Neumann term originally evaluated `stokesEq:velocity()`
   on Membrane elements, outside the Stokes exporter domain.  It is now
   assembled from the Lumen side into the shared Apical Darcy trace DoFs.
6. Artificial Stokes velocity/pressure extensions into porous subsets were
   removed.  r1 fell from 54,366 to 36,646 DoFs.
7. Refined Uzawa GMG and full ILU produce NaNs.  Element-patch GMG is finite
   but diverges/stagnates.  r1/r2 therefore use pivoted SuperLU plus one
   residual-correction solve; this is an algebraic diagnostic, not the final
   scalable solver.

## Reproduction

```text
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell -ex \
  code/flow/qi_2025_paper_interface_stokes_darcy_simple_long_tube.lua \
  -numRefs N -output /tmp/qi_compact_rN
```

The exact r1/r2 trial prefixes used during this checkpoint were
`/tmp/qi_compact_direct_r1` and `/tmp/qi_compact_direct_r2`.

## Algebraic and physical results

| Level | Elements | DoFs | Max component defect | Lumen relative balance | Darcy exit mismatch |
|---:|---:|---:|---:|---:|---:|
| r0 | 11,045 | 5,012 | 5.37e-14 | 4.46e-13 | 97.01% |
| r1 | 88,360 | 36,646 | 7.10e-14 | 7.17e-13 | 79.56% |
| r2 | 706,880 | 278,763 | 5.15e-16 | 3.44e-15 | 99.29% |

The Darcy exit mismatch compares the Apical flux imposed by the Stokes field
with the OuterWall flux reconstructed from the Darcy pressure gradient.  The
1% acceptance criterion fails on every level and does not improve
monotonically.

The lumen split is also non-convergent.  At r0 the Outlet flux is
`-1.733322e-3` and Apical flux is `2.214379e-3` mm^3/s; at r1 they are
`-2.186360e-6` and `4.775460e-4`; at r2 they are `9.264574e-9` and
`4.723720e-4` mm^3/s.  Per-level lumen balance remains excellent but does not
establish coupled physical correctness.

## Exported membrane/inter field trend

| Level | Membrane mean speed | Membrane max speed | Inter mean speed | Inter max speed |
|---:|---:|---:|---:|---:|
| r0 | 2.1749e-3 | 1.1611e-2 | 3.6950e-4 | 5.4819e-3 |
| r1 | 5.8127e-3 | 2.2425e-1 | 4.6655e-4 | 2.6366e-2 |
| r2 | 8.5998e-5 | 3.4676e-3 | 1.4511e-5 | 1.2595e-3 |

Speeds are in mm/s.  The non-monotone changes rule out mesh convergence.
Darcy pressure ranges likewise spike at r1 and collapse at r2.

## ParaView inspection

At the user's explicit request, the following unverified r2 debugging datasets
were opened in ParaView 5.11.2:

- `/tmp/qi_compact_direct_r2_lumen_t0000.vtu`
- `/tmp/qi_compact_direct_r2_membrane_t0000.vtu`
- `/tmp/qi_compact_direct_r2_inter_t0000.vtu`

They must not be presented as verified results or copied to `Results_final`.

## Next technical target

Audit the Darcy Neumann/interface row against the FV1 conservative boundary
flux, preferably with a minimal manufactured straight-interface test that
exposes the assembled Darcy boundary flux directly.  Do not tune permeability,
BJS coefficients, or physiological pressures to hide the failure.  After the
interface row is conservative, restore a scalable block/GMG solver and repeat
r0--r2 flow, pressure-drop, and velocity-profile convergence checks before
moving to nephron geometry.
