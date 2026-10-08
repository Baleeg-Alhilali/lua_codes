# Anisotropic projector experiments — 2026-09-29

## Current retained configuration

- Active build: `/Users/alhilaba/UG4_promesh_ogrid`
- The experimental radial corrections described below were removed.
- The retained projector is the previously verified progressive projector.
- Its centreline/surface displacement limit is `0.002 * parentScale` (0.2%).
- Default mesh-generation anisotropy is now `1`.
- `-anisotropy` can still explicitly override the default.

## Observed problem

With a coarse mesh generated using anisotropy 4, a vertex inserted halfway
between two axial rings is initially the Cartesian midpoint of its parent
edge. On a curved tube this midpoint is closer to the interpolated centreline
than the two parent vertices. The new intermediate ring therefore has a
slightly smaller radius. This contributes to visible sharp/constricted rings.

The old straight-edge projector and the 0.2% progressive projector look almost
identical because the 0.2% cap also limits the outward radial correction.

## Controlled comparison meshes

All comparison runs used the third trace, radius 1 um, membrane thickness
0.2 um, O-grid 16, BoxRefs 0, and two refinements unless stated otherwise.

- `runs/third_trace_anisotropy4_before_axial_projector_fix*.ugx`
  - Legacy straight-parent-edge axial behavior.
  - Refinement 1: 0 degenerate volumes.
  - Refinement 2: 0 degenerate volumes.
- `runs/third_trace_anisotropy4_after_progressive_projector_fix*.ugx`
  - Retained 0.2% progressive behavior.
  - Refinement 1: 0 degenerate volumes.
  - Refinement 2: 0 degenerate volumes.

## Radial-correction approach tested

For a true axial edge, the experiment calculated:

1. the Cartesian midpoint of the two surface/radial-layer vertices;
2. the midpoint of their two exact SWC centreline positions;
3. the contracted midpoint radius;
4. the prescribed exact radius at the averaged axial parameter;
5. a corrected radius
   `linearRadius + factor * (exactRadius - linearRadius)`;
6. a separately capped centreline displacement.

This separates radial contraction from centreline sagitta, but moving vertices
after tetrahedral refinement can still create poor or degenerate membrane
tetrahedra.

## Test results

- Full radial correction applied broadly to axial parents:
  - 14 degenerate membrane volumes at refinement 1.
- Full correction restricted to true axial edges:
  - 8 degenerate membrane volumes at refinement 1.
- 50% radial correction:
  - refinement 1 passed;
  - refinement 2 produced 1,429 degenerate membrane volumes when repeated.
- 50% correction only when no progressively corrected parent was detected:
  - refinement 1 passed;
  - refinement 2 produced 1 degenerate membrane volume.
- 40% correction with the same staging:
  - refinement 1 passed;
  - refinement 2 produced 6 degenerate membrane volumes.
- Separate persistent `radius corrected` marker with 50% correction:
  - refinement 1 passed;
  - refinement 2 produced 65 degenerate membrane volumes.
- 25% staged correction:
  - refinement 1: 0 degenerate volumes;
  - refinement 2: 0 degenerate volumes;
  - files: `runs/third_trace_anisotropy4_first_ring_25pct_radial*.ugx`.
- 25% staged correction with membrane thickness 0.5 um:
  - refinement 1: 0 degenerate volumes;
  - refinement 2: 0 degenerate volumes;
  - files: `runs/third_trace_anisotropy4_membrane_0p5_radius_correction*.ugx`.

Although the 25% cases passed orientation validation, their visual improvement
was insufficient. Therefore this experimental correction was not retained.

## Recommended future direction

Do not increase the post-refinement projection displacement globally. A 0.5%
full projection test performed earlier produced 324 degenerate membrane
volumes at refinement 2.

The robust solution should create intermediate axial rings parametrically from
the SWC centreline, transported frame, and prescribed radial coordinate before
tetrahedralization, then construct/refine the volumes around those rings. That
keeps every ring radius consistent and moves the interior coherently, instead
of moving already-connected tetrahedral vertices independently afterward.

An alternative is a quality-aware deformation/backtracking stage that projects
the boundary and moves the interior elastically while rejecting any step that
reduces a tetrahedron Jacobian below a safe threshold.

