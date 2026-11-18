# Nephron mesh geometry pipeline: tests, fixes, and current baseline

Consolidated: **2026-10-08**  
Purpose: preserve what was tried, what failed, what was retained, and the exact
baseline from which the next geometry improvements should start.

## 1. Current recommended baseline

The geometry baseline is the dated snapshot:

`/Users/alhilaba/Desktop/O-grid_W_surface_fixed_2026-10-04`

It is separate from, and did not overwrite:

`/Users/alhilaba/Desktop/O-grid_working_FINAL`

The current verified production choices are:

| Setting | Production value | Reason |
|---|---:|---|
| Anisotropy | `1` | Values above 1 produced poor/distorted refined elements and unequal-looking intermediate rings. |
| Circumferential O-grid vertices | `10` | Reduced volume/vertex count while retaining acceptable cross-sections. |
| BoxRefs | `4` | Good OuterWall resolution without the cost of BoxRefs 5. |
| Lumen coarse radial topology | center + midpoint ring + Apical ring | Removing the midpoint ring created alternating coarse-vertex diagonals during refinement. |
| Refinement used for production/GMG | global conforming, `crossOnly=0` | Lumen, Membrane, Inter, and OuterWall remain conforming and retain a usable hierarchy. |
| SWC resampling | `0.5 um` | Stable spatial sampling used in the verified build. |
| Pre-smoothing | `10 degrees`, strength `0.85`, `50` passes | Removed the W-shaped refined surface while allowing exact projected placement. |
| Projector frame transport | deterministic rotation-minimizing frame, `1024` steps | Prevents frame twisting and makes results reproducible. |
| Lumen radius | `1.0 um` | Current test/production radius. |
| Membrane thickness | `0.2 um` | `0.5 um` was tested but did not correct the anisotropy/refinement problem. |
| OuterWall port transition | BoxRefs-cell-centered, triangular, aligned adaptive transition | Avoids long edges and avoids a separate artificial small square. |

Final geometry validation for this baseline:

- coarse mesh: **29,211 vertices**, **73,591 volumes**;
- refinement level 1: **0 indeterminate/degenerate volumes**;
- refinement level 2: **0 indeterminate/degenerate volumes**;
- maximum Apical/Basolateral radial residual at levels 1 and 2:
  **4.15097e-10 mm**;
- anisotropy 1, O-grid 10, BoxRefs 4, global conforming refinement.

Primary inspection meshes:

- `work/ring_smooth10_l2.ugx`
- `work/ring_smooth10_l2_refined_1.ugx`
- `work/ring_smooth10_l2_refined_2.ugx`
- `runs/3rd_trace_new_build_test.ugx`
- `runs/3rd_trace_new_build_test_refined_1.ugx`
- `runs/3rd_trace_new_build_test_refined_2.ugx`

## 2. Original Neuro Collection pipeline

The original pipeline imported an SWC centerline and built a structured tube
inside a padded box:

1. Read the SWC samples and radii.
2. Smooth/resample the centerline.
3. Build oriented polygonal rings along the centerline.
4. Connect corresponding ring vertices axially.
5. Create Lumen and Membrane volumes from the nested rings.
6. Create a surrounding box and tetrahedralize the Inter region with TetGen.
7. Assign Lumen, Membrane, Apical, Basolateral, Inlet, Outlet, Inter, and
   OuterWall subsets.
8. Serialize `npSurfParams` and NeuriteProjectors into the coarse UGX.
9. Use the projectors when new surface vertices are created during refinement.

The first improved O-grid contained:

- one centerline vertex per station;
- an internal Lumen ring;
- the Apical ring;
- an internal Membrane ring;
- the Basolateral ring.

This gave a prismatic core plus hexahedral outer Lumen and Membrane layers. It
also ensured that both volumes had internal vertices rather than only boundary
vertices.

## 3. Early refinement fixes retained in `O-grid_working_FINAL`

The following issues were fixed before the OuterWall-port work:

- terminal caps were detected using actual Lumen ownership and axial position;
- internal Lumen faces were excluded from Inlet/Outlet;
- Apical rim edges kept Apical priority when also touching a cap;
- all other terminal cap edges became Inlet or Outlet;
- `npSurfParams` were propagated through Inlet/Outlet refinement;
- orphan child objects were removed between generated levels;
- element orientation was checked and repaired after each refinement;
- detached objects no longer became parents of the next level;
- flat `_refined_N.ugx` inspection files were written instead of visually
  stacking all parent levels.

Earlier verified O-grid test:

- input: `3rd_trace_243_third.swc`;
- O-grid 16, anisotropy 2, BoxRefs 4;
- level 1: zero degenerate volumes;
- level 2: 8 ordinary orientation corrections and zero degenerate volumes;
- the former 96 collapsed terminal tetrahedra were eliminated.

This build was a useful reference, but anisotropy 2 was later rejected for the
new OuterWall-port geometry.

## 4. Anisotropy and membrane-thickness investigation

### Observed defect

With anisotropy 4, the refined intermediate axial rings appeared slightly
smaller than the coarse rings. The geometry before/after refinement was close,
but the unequal radius produced visible surface waviness and poor elements.

### Tests

| Test | Result | Decision |
|---|---|---|
| anisotropy 4 | Intermediate rings did not match coarse-ring radius; bad elements appeared. | Rejected. |
| other anisotropy values above 1 | Similar refinement-quality problems. | Rejected. |
| membrane thickness `0.5 um` | Did not bring the geometry closer to the desired result. | Rejected as a fix. |
| anisotropy 1 | Most stable element geometry. | Retained and made the working default/required value. |

The experiments were preserved for later study, but production returned to
anisotropy 1.

## 5. Extending Inlet and Outlet to the OuterWall

### Required behavior

- Project both terminal sections onto their nearest OuterWall face.
- Extend the SWC before mesh generation rather than stretching only the mesh.
- Make the final cap parallel to, and share the normal of, the selected wall.
- Use the least-curved path and do not continue in the direction of the rest of
  the nephron when that produces unnecessary bending.

### Implemented behavior

- The smoothed, millimetre-scaled SWC is extended at both ends.
- A low-bending cubic transition is followed by a straight wall-normal tail.
- Fixed OuterWall bounds are stored in SWC metadata so the extension does not
  move the box again.
- The final cross-section is snapped exactly to the chosen box plane.
- The terminal Basolateral rim becomes the exact hole in the OuterWall; no wall
  face is left over the opening.
- Inlet/Outlet/TripleJunction vertices are snapped back to the wall plane after
  each projected refinement.

### Curvature/padding test

- Requested 10 um padding gave root/tip maximum curvatures of approximately
  391/mm and 271/mm. TetGen crashed on the folded/intersecting boundary.
- 30 um padding reduced them to approximately 123/mm and 85/mm and meshed.
- A minimum effective terminal padding of 25 outer radii was retained through
  `terminalPaddingFactor=25`.

The unrelated value `0.06` in the extension search is a candidate step used to
find a low-curvature route. It is not a 6% NeuriteProjector displacement.

## 6. OuterWall and portal-refinement experiments

### Surface element choice

- OuterWall quads/squares caused awkward portal transitions.
- A triangular OuterWall was retained.
- Prisms were considered, but a prism is a 3D volume and cannot replace a 2D
  wall face directly. A prism boundary layer remains a separate future option,
  not a fix for the observed surface problem.

### Iteration history

1. **Refine the complete portal-bearing box face**
   - Removed some unrefined interiors but introduced unnecessary global cost and
     long connections.
   - Rejected.

2. **Small refined square around the circular port**
   - Improved local density but the small square itself created difficult
     transitions and corner-biased edges.
   - Rejected as the final topology.

3. **Local 2x2/3x3 neighborhood**
   - Refined only nearby BoxRefs cells and reduced long rim-to-wall edges.
   - Useful intermediate step.

4. **Adaptive 4-to-2-to-1 transition**
   - Port cell edges received four segments, neighbors two, distant cells one.
   - Shared edge vertices and cell-center triangulation made the transition
     conforming.
   - Valid test: 0 open PLC edges, 0 non-manifold edges, and 0 degenerate
     volumes after one refinement.

5. **No-small-square variant**
   - Removed the artificial square and its clearance enlargement.
   - Refined the actual BoxRefs cell containing the circular port.
   - Coarse OuterWall: 496 vertices, 1,500 edges, 992 triangular faces;
     one refinement had zero degenerate volumes.

6. **Cell-centered terminal variant**
   - Moved the extension target to the exact center of the selected structured
     wall cell.
   - Removed experimental cross constraints and avoided creating a second box.
   - Worked with BoxRefs 5; BoxRefs 3 caused a deterministic TetGen crash.

7. **Reduced aligned 2-to-1 transition**
   - BoxRefs changed to 4; port-cell boundary reduced to two segments.
   - Square-boundary vertices were connected to angularly aligned rim vertices.
   - Coarse: 51,095 vertices, 168,359 volumes, 0 orientation corrections;
     ref1: 12 corrections, 0 degenerate volumes.

8. **Final aligned 4-to-2-to-1 variant**
   - Kept BoxRefs 4.
   - Restored four segments in the port cell and kept explicit angular edge
     assignment.
   - Coarse: 51,262 vertices, 169,249 volumes, 0 orientation corrections;
     ref1: 15 corrections, 0 degenerate volumes.

The final requested visual adjustment retained BoxRefs 4 and the improved edge
assignment.

## 7. Reducing the coarse nephron mesh

### O-grid count

The circumferential count was reduced to 10 to make a one-level refined mesh
comparable in cost to two or three refinements of the older coarse mesh.

### Removed Lumen midpoint ring experiment

`coarseLumenCenterOnly=1` removed the ring between the center and Apical ring.
It reduced vertices, but each full center-to-Apical sector then had to be
triangulated. Refinement necessarily created diagonals toward alternating old
coarse angular vertices, producing the unwanted faces visible at the inlet and
through the Lumen/Membrane.

Decision:

- restore the midpoint support ring (`coarseLumenCenterOnly=0`);
- keep the center-only form only as an explicit diagnostic option.

This was a topology issue, not a projector-position issue.

## 8. Cross-sectional-only refinement investigation

### Goal

Refine the nephron cross-section with the NeuriteProjector without creating new
axial stations, while refining Inter and OuterWall normally.

### What worked

- Projector vertices were created at reasonable surface positions.
- Cross-sectional masks prevented new axial sections.
- A local topology test showed that midpoint-to-midpoint connections could be
  produced after explicit face/volume mask handling.

### Problems found

- New vertices sometimes connected diagonally to an old coarse-ring vertex
  rather than to the matching midpoint on the next ring.
- Basolateral angular midpoints at the terminal could remain on a polygon chord
  with zero neurite parameters.
- Zero parameters propagated into the next refinement and collapsed cells.
- Partial refinement produced constrained faces that NavierStokesFV1 treated as
  a non-regular grid.
- The attractive `*_active_level` display was only the active refined surface;
  it did not mean the entire mixed nephron/Inter/OuterWall grid was conforming.

### Cross-only connectivity fixes implemented/tested

- Structured Lumen/Membrane hexahedra were preserved.
- Two TetGen Basolateral triangles were represented as constrained children of
  the structured Basolateral quad.
- Matching face and volume masks were marked so adjacent rings connected
  midpoint-to-midpoint.
- Node-dependency order 1 was retained.
- An explicit nonzero local mask was preserved when constraint propagation
  promoted a general refinement mark.
- Missing neurite parameters were restored on coincident constrained copies.

Verified diagnostic configuration:

`anisotropy=1, ogrid=10, BoxRefs=4, crossOnly=1, numRefs=2`

- no new axial sections at levels 1 or 2;
- zero degenerate volumes at both levels;
- topology check found no angular midpoint-to-old-ring diagonal.

### Production decision

Cross-sectional-only refinement remains experimental. It was not selected for
the coupled conforming mesh because Lumen/Membrane partial refinement and
regular Inter/OuterWall refinement did not give the simple globally conforming
hierarchy required by the multigrid workflow.

Production returned to `crossOnly=0`: refine everything globally and
conformingly, then revisit cross-sectional refinement later.

## 9. Boundary subset corrections

### Basolateral persistence

Refined longitudinal outer child faces are recognized by their complete outer
radial neurite parameters and reassigned to Basolateral. This keeps the
Membrane--Inter interface present and projected at every level.

### MembraneWall

The terminal membrane annulus lying on the box plane receives its own
`MembraneWall` subset. It is not the internal longitudinal Basolateral boundary,
because it will require a different boundary condition.

`MembraneWall` uses the neurite projector during child creation and is snapped
back to the wall plane afterward.

### Terminal repair

Older coarse meshes mislabeled one-sided terminal annuli as Basolateral. A
repair pass reclassifies those terminal faces/edges before refinement. This
prevents unsupported Basolateral DoFs and incorrect transport coupling.

## 10. Hierarchy and file-format correction

An earlier attempt serialized the visible hierarchy directly into refined UGX
files and created `*_active_level` companions. That was not the workflow used
by the known working multigrid build.

Correct behavior:

- the coarse UGX stores the projector and `npSurfParams`;
- generator-side `_refined_N.ugx` files are standalone flat inspection meshes;
- the solver loads the coarse UGX;
- `RefineAndRebalanceDomain` creates the real parent/child hierarchy in memory;
- the approximation space initializes all in-memory levels for GMG.

Therefore, flat refined files are for geometry inspection and validation. They
must not replace the projector-backed coarse mesh as the GMG input.

## 11. Projector-radius and W-shaped surface investigation

### Reference audit

The September reference `4th_trace_third.ugx` was inspected directly:

- 91,233 Inter tetrahedra;
- 9,248 Membrane hexahedra;
- Lumen: 4,624 hexahedra plus 4,624 prisms;
- no pyramid volumes;
- each refined volume had eight children;
- axial child vertices used Cartesian parent midpoints;
- transverse children were projected.

The old mesh appeared smooth because its centerline bends were milder, not
because it used a different exact-radius axial projector.

### Projection-cap experiments

| Experiment | Outcome |
|---|---|
| Full exact independent surface projection before coarse-centerline repair | 18 invalid volumes in two high-curvature intervals. |
| 15% move cap | 22 degenerate Membrane volumes. |
| 5% move cap | 2 degenerate Membrane volumes. |
| 4% move cap | 2 degenerate Membrane volumes. |
| 3% move cap | Level 1 valid, but 30,820/54,860 vertices limited; max residual `0.000186847`. |
| Exact surface parents, limited interior parents | 71 degenerate volumes. |
| Coherent edge/face/volume exact map | 18 degenerate volumes near `t=0.821--0.825`. |
| Chord-centered normalization | Reduced invalid count only to 15. |
| Curvature-local limiter | Level 1 valid, residual `0.0000889515`; level 2 produced 292 invalid volumes. |
| Inherited local limiter flags | Still 234--402 invalid level-2 volumes. |
| Curvature-aware coarse subdivision | Added one station, then TetGen crashed. |

The invalid exact-projection cells occurred where adjacent 0.5 um samples
turned by approximately 22.7 and 25.5 degrees. The midpoint rule merely hid the
coarse-centerline defect by pulling new rings inward.

### Final retained repair

1. Smooth sharp control-point turns above 10 degrees for 50 passes.
2. Use the deterministic 1024-step rotation-minimizing frame.
3. Apply exact projector positions to both axial and transverse children.

This removed the W-shaped Apical/Basolateral surface and achieved numerical
radius accuracy while keeping all volumes valid through two global levels.

## 12. Flow-mesh and solver validation

### Whole geometry

Generated files:

- `runs/3rd_trace_new_build_test.ugx` (coarse whole domain);
- `runs/3rd_trace_new_build_test_refined_1.ugx`;
- `runs/3rd_trace_new_build_test_refined_2.ugx`.

Settings: anisotropy 1, O-grid 10, BoxRefs 4, global refinement, pre-smoothing
10 degrees/50 passes, radius 1 um, membrane 0.2 um.

### Why a Lumen companion was needed

The whole mesh contains constrained Membrane--Inter interface objects.
NavierStokesFV1 reports `only regular grid implemented` even if equations are
defined only on Lumen. No NavierStokes plugin change was retained.

`extract_lumen_new_build.lua` removes Membrane and Inter volumes and unused
constrained objects, while retaining Lumen, Apical, Inlet, Outlet,
TripleJunction, all serialized NeuriteProjectors, and `npSurfParams`.

Result:

`runs/3rd_trace_new_build_test_lumen.ugx`

### Solver tests and fixes

| Attempt | Result |
|---|---|
| Scalar Gauss--Seidel GMG | Diverged. |
| ComponentGaussSeidel pressure grouping | Failed in `extract_by_grouping`. |
| Element Vanka with iterative coarse solve | Fine smoothing worked; coarse solve failed. |
| Gathered direct solve with the wrong transient/RAP setup | Impractically slow. |
| Temporary force-regular bridge experiment | Reverted; no retained bridge/core change from this solver test. |
| Proven old-flow configuration | Successful. |

The successful configuration copied the earlier working flow logs:

- stationary Stokes FV1;
- damped ILU (`0.7`) on refined levels;
- rediscretized multigrid operators (`RAP=false`);
- `StdTransfer`;
- gathered SuperLU only on level 0;
- all requested projector levels created before MPI distribution.

Verified results:

| Run | Top-level DoFs | Linear iterations | Solve time | Result |
|---|---:|---:|---:|---|
| serial ref1 | 321,276 | 6 | 5.926 s | PASS |
| serial ref2 | 2,360,700 | 8 | 67.106 s | PASS |
| MPI 4 ref1 | 321,204 | 13 | 2.851 s | PASS |
| MPI 4 ref2 | 2,360,484 | 13 | 31.202 s | PASS |
| MPI 8 ref1 | 321,204 | 13 | 2.084 s | PASS |

The flow script now writes `velocity`, `u`, `v`, `w`, and `pressure`. Output
names begin with `new_build_` and target:

`/Users/alhilaba/Documents/Promesh_swc/Model/Results`

## 13. Current pipeline, end to end

1. Select the original SWC.
2. Pre-smooth high-curvature control points using the verified 10-degree,
   strength-0.85, 50-pass settings.
3. Resample at 0.5 um.
4. Override the SWC radius to outer radius = lumen radius + membrane thickness.
5. Extend both ends to selected OuterWall cell centers with low-curvature cubic
   transitions and wall-normal tails.
6. Preserve fixed box bounds in metadata.
7. Build 10-vertex O-grid stations with a deterministic rotation-minimizing
   frame.
8. Keep the Lumen midpoint support ring.
9. Build structured Lumen and Membrane volumes.
10. Build the BoxRefs-4 triangular OuterWall with locally aligned adaptive port
    transitions.
11. Use terminal Basolateral rims as portal holes.
12. Tetrahedralize only Inter and validate the PLC first.
13. Remove overlaps, duplicates, detached objects, and repair orientations.
14. Assign biological subsets, including `MembraneWall` at terminal annuli.
15. Serialize `npSurfParams` and projectors in the coarse UGX.
16. For geometry inspection, generate flat global refined levels and validate
    radius/orientation at every level.
17. For multigrid, ignore flat refined files: load the coarse mesh and generate
    the hierarchy in memory.
18. For Lumen-only flow, derive the regular lumen companion in Lua and solve
    stationary Stokes with the validated GMG configuration.

## 14. Known limitations and future work

- Cross-sectional-only refinement is not yet the production coupled-mesh mode.
  It needs a conforming transition strategy compatible with Inter/OuterWall and
  GMG, without constrained objects visible to regular-grid-only discretizations.
- A prism boundary layer in Inter is possible but is a new 3D meshing feature,
  not a replacement for triangular wall faces.
- Very coarse BoxRefs values can make portal PLC reconstruction fail.
- Changes to smoothing/resampling must be tested at both refinement levels;
  level-1 success alone was not predictive.
- Independent per-vertex projection should never be increased without checking
  coherent cross-section deformation and signed volume at levels 1 and 2.
- Preserve anisotropy 1 unless a new method passes geometry and solver tests.
- Do not remove the Lumen midpoint support ring merely to reduce the mesh count.
- Do not edit NavierStokes or UG4 core for new solver work. The later workflow
  constraint is to make new changes only in neuro_collection and Lua. Earlier
  geometry snapshots contain historical projector/hanging-node changes and must
  be treated as fixed dependencies, not reopened casually.

## 15. Validation checklist for every future change

1. Generate a new uniquely named coarse mesh; never overwrite a dated snapshot.
2. Check PLC open and non-manifold edges before TetGen.
3. Check coarse orientations and zero-volume cells.
4. Generate global refinement level 1 from the coarse mesh.
5. Check subset counts, projector parameters, orientation, signed volume, and
   surface-radius residual.
6. Repeat at level 2; do not accept a change based only on level 1.
7. Inspect Inlet, Outlet, Apical, Basolateral, MembraneWall, and OuterWall
   separately in ProMesh.
8. Confirm midpoint-to-midpoint topology across adjacent rings.
9. Confirm no coarse-parent overlay is being mistaken for active refined faces.
10. Build the in-memory hierarchy from the coarse UGX and run stationary Stokes
    at ref1 and ref2.
11. Run serial and MPI tests and require `FLOW_TEST_SUCCESS`.
12. Record the command, source hashes, mesh counts, residuals, solver iterations,
    and output paths in this log or a dated follow-up.

## 16. Reference files and snapshots

- `work/BUILD_SNAPSHOT_2026-10-04.md`
- `work/outerwall_portal/README.md`
- `work/outerwall_portal/REFINEMENT_CONNECTIVITY_FIX_2026-10-04.md`
- `work/FLOW_NEW_BUILD_VALIDATION_2026-10-05.md`
- `work/SCALE_ADAPTIVE_GEOMETRY_2026-10-08.md`
- `work/test_ring_projection_audit.lua`
- `work/ring_smooth10_l2.log`
- `/Users/alhilaba/Desktop/O-grid_working_FINAL`
- `/Users/alhilaba/Desktop/O-grid_W_surface_fixed_2026-10-04`

Current key SHA-256 values at consolidation time:

```text
cd1e061beed84ac3db9c2fb418428b779fd09e1a5f5c800b5b6cb66e9bc8fa8f  active neurites_from_swc.cpp
613d6054ea7d80c18ce7252cefcb7e684fd91723d5be2be2dc87850aec327315  active neurites_from_swc.h
291141f46c17936874dec1061f7752c60de91aecdbccdf94a8d73e4eec5d463b  active neuro_collection_plugin.cpp
204cbce106117fd44ec1f577376863f7e3fb8cefda09e765fa64dc39b3cb2c7e  stable neurite_projector.cpp
b4fb1b1a7f3cf9d43eafd5a3a8397f66f4ab7109eb15d6d541a46d9f6875cb7f  runs/3rd_trace_new_build_test.ugx
bed27b29cff31559488aa22dfdcf06b8e8383dce7c543c70864433fd56f11487  runs/3rd_trace_new_build_test_lumen.ugx
877cb3f0841a4baf1246ec300e5a49052eca052578b9f4ed12a8bd5419dd84ef  flow_projector_new_build.lua
```

## 17. Variable-radius and thick-membrane audit (2026-10-08)

The geometry pipeline was audited for scale-dependent thresholds. Numerical
position tolerances in the Neuro Collection additions now derive from the
local radius, element scale, grid diagonal, or OuterWall diagonal. Terminal
sampling and box padding now derive from the requested outer radius.

The audit also separated a numerical-threshold failure from a physical
geometry failure: the third trace has a minimum discrete centerline curvature
radius of about 5.115 um, so an unscaled tube with 20--25 um outer radius
self-intersects. A new preflight rejects this before TetGen and reports the
required `-centerlineScale`.

Validated outputs:

- reference regression: `runs/scale_regression_r1_t0p2.ugx`
- 10 um lumen radius, 15 um membrane, centerline scale 6:
  `runs/scale_valid_r10_t15_s6_final.ugx`
- both generated one projected refinement with zero indeterminate/degenerate
  volumes.

Full parameter definitions, threshold classification, commands, counts, and
source hashes are in `work/SCALE_ADAPTIVE_GEOMETRY_2026-10-08.md`.

Axis-specific OuterWall clearance was subsequently added through `-paddingX`,
`-paddingY`, and `-paddingZ`. The controls were validated with 40/50/70 um
X/Y/Z clearances in `runs/axis_padding_test.ugx`.

## 18. Adaptive Inter and bend-safe anisotropy (2026-10-08)

The coarse Inter now receives a scale-relative distance field: tetrahedra are
small near Basolateral, grow through a transition band, and remain coarse near
OuterWall. Coarse nephron station spacing now uses anisotropy 4 on straight
paths and recursively approaches anisotropy 1 only where spline/chord error
detects a bend. The accepted large-radius mesh has 37,145 coarse volumes; one
global projected refinement completed with zero indeterminate/degenerate
volumes. The detailed algorithm, rejected over-refined trial, parameters,
measured edge grading, two-level global-refinement validation, and rollback flags are recorded in
`work/ADAPTIVE_INTER_AND_ANISOTROPY_2026-10-08.md`.
