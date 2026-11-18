# Cross-sectional projector connectivity fix — 2026-10-04

## Defect

Cross-sectional refinement created the projector vertices at plausible positions,
but the refined membrane topology joined some new vertices diagonally to an old
coarse-ring vertex.  At the terminal OuterWall, the Basolateral angular midpoint
was also left on the polygon chord and acquired zero neurite surface parameters.
On a second refinement those zeros were averaged into the membrane and collapsed
some terminal hexahedra toward the centerline.

## Implemented fix

- Preserve the structured Lumen/Membrane hexahedra.
- Represent TetGen's two exterior Basolateral triangles as constrained children
  of the original structured Basolateral quadrilateral.
- Mark matching local face and volume masks so transverse midpoints on adjacent
  rings connect midpoint-to-midpoint; axial edges are not split.
- Keep node-dependency order 1 enabled for the structured local refinement.
- In the hanging-node refiner, retain an explicit nonzero local mask when
  constraint propagation promotes the general refinement mark.
- Classify the exact terminal outer membrane rim as Basolateral while leaving
  the surrounding planar box faces as OuterWall.  Consequently only the circular
  rim receives the NeuriteProjector.
- Retain a safety pass that restores neurite parameters on geometrically
  coincident constrained-interface copies, if such copies occur.

## Verified configuration

`-anisotropy 1 -ogrid 10 -BoxRefs 4 -crossOnly 1 -numRefs 2`

Results:

- refinement levels 1 and 2 introduced zero new axial sections;
- level 1: zero degenerate/indeterminate volumes;
- level 2: zero degenerate/indeterminate volumes;
- terminal outer-rim vertices remain in Basolateral, not OuterWall;
- the one-level topology check found straight matching midpoint connections and
  no angular-diagonal midpoint-to-coarse-ring connections.

Saved source snapshots in this directory:

- `neurites_from_swc.cpp`
- `neurite_projector.cpp`
- `hanging_node_refiner_base.cpp`

Inspection meshes are in:

`/Users/alhilaba/Documents/Promesh_swc/Model/runs/third_trace_projector_edges_fixed*.ugx`
# Basolateral children, wall contact, and flat refinement output

The terminal membrane annulus that lies on an OuterWall plane now uses a
dedicated `MembraneWall` subset.  It is not part of the longitudinal
`Basolateral` transport boundary, so the two locations can receive different
boundary conditions.  `MembraneWall` uses the neurite projector while children
are created and is then snapped exactly back to its OuterWall plane.

Cross-sectional hanging-node refinement initially creates the longitudinal
outer child faces as sides of `Membrane` volumes.  A post-refinement pass now
recognizes active faces whose complete vertex set has the outer radial neurite
parameter and restores those faces, their edges, and their vertices to
`Basolateral`.  This makes the refined Basolateral interface persistent at every
level and ensures that its next children use the same projector.

The earlier hierarchy serialization was incorrect for this workflow.  As in
`O-grid_working_FINAL`, every `_refined_N.ugx` is now a standalone flat mesh
containing level N only.  No `_active_level` companion is written.  The flow
solver loads the projector-backed coarse UGX and creates the required hierarchy
in memory with `RefineAndRebalanceDomain`; it does not load the pre-generated
flat `_refined_N.ugx` file for geometric multigrid.

The lumen-midpoint support ring is enabled again by default
(`-coarseLumenCenterOnly 0`).  The center-only coarse lumen removed this ring
and forced refinement to triangulate each full center-to-Apical sector.  Those
children necessarily retained diagonals toward alternating coarse angular
vertices, producing the visible pattern at the inlet and through the lumen.
Keeping one midpoint ring restores the original nested radial refinement
connectivity without changing the current OuterWall-port construction.  The
smaller center-only experiment remains selectable with
`-coarseLumenCenterOnly 1`.

## Production refinement mode for multigrid

The production default is global conforming refinement (`-crossOnly 0`).  The
Lumen, Membrane, Inter, and OuterWall are refined together, so every interface
is conforming.  Generator-side `-numRefs` writes flat meshes for inspection and
validation.  Solver-side `-numRefs` is separate: `Dif_itterative_flow_projector.lua`
loads the coarse mesh, calls `RefineAndRebalanceDomain`, and thereby builds the
true geometric-multigrid hierarchy in memory before initializing the
approximation-space levels.

Cross-sectional-only refinement remains available explicitly with
`-crossOnly 1` for future development.

## Axial projection audit

The projector was audited directly against the saved
`~/Desktop/O-grid_working_FINAL` source and executable.  That verified
reference does not use the later `0.002 * parentScale` progressive
displacement.  It projects transverse/circumferential children exactly, while
a child whose parent spans two axial stations remains at the Cartesian parent
interpolation.  Descendants inherit this state with
`NP_AXIAL_INTERPOLATION_FLAG`, so every level treats the complete axial ring
consistently.  This reference rule is restored in the active build.  The
unrelated `0.06` value in the SWC extension script is only a candidate-step
increment for choosing a low-curvature path to the wall.

Removing the cap and using the exact projected node was tested with anisotropy
1, O-grid 10, BoxRefs 4, and one global refinement.  It produced 18 degenerate
volumes (14 Membrane and 4 Lumen), so exact independent vertex projection was
not retained.  The restored reference projector passes one global refinement
with zero degenerate volumes, both before and after extending the terminals to
OuterWall.  Making newly inserted axial rings converge to the exact analytical
radius beyond this reference behavior still requires coherent cross-section
deformation or parametric ring construction, not a larger independent
per-vertex move.

## Radius/hierarchy experiments after the audit

The reported Apical/Basolateral hills and valleys are real geometry, not an
inactive-parent display artifact.  A numerical audit now compares every new
surface vertex with the serialized NeuriteProjector and reports the radial
residual before orientation validation.

The following alternatives were built and tested on the terminal-to-OuterWall
third trace (`anisotropy=1`, O-grid 10, BoxRefs 4, global refinement):

- A global 15% projection-move cap created 22 degenerate Membrane volumes.
- Caps of 5% and 4% each left 2 degenerate Membrane volumes.
- A global 3% cap passed level 1, but still limited 30,820 of 54,860 checked
  surface vertices; maximum exact-ring residual was `0.000186847`.
- Exact projection of constant-radius surface parents while limiting interior
  parents created 71 degenerate volumes.
- One coherent exact map for edge, face, and volume children reduced this to 18
  degenerate volumes, all in two tight axial intervals near 0.821--0.825.
- Chord-centered radius normalization reduced the count only to 15; it did not
  make those sharp-bend cells valid.
- A curvature-local level-1 limiter passed level 1 with zero degenerates and
  maximum residual `0.0000889515`, but level 2 created 292 degenerate volumes.
  Inheriting the local flag still left 234--402 degenerate level-2 volumes,
  because other level-1 spline-projected cells were positive but non-convex.
- Enabling the existing curvature-aware coarse subdivision inserted only one
  extra axial segment (610 to 611), but TetGen crashed while tetrahedralizing
  Inter.  This experimental path remains disabled.

The active production rule therefore matches the stable
`O-grid_working_FINAL` hierarchy: every axial child uses the Cartesian parent
interpolation and carries `NP_AXIAL_SAFETY_FLAG` to descendants; transverse
children still use the exact projector.  It passed two global levels:

- level 1: 0 degenerate volumes (15 orientation corrections);
- level 2: 0 degenerate volumes (13 orientation corrections).

For the tested 1 micrometre lumen radius, the safe rule's level-1 radial audit
had mean residual `0.000019058` and maximum residual `0.000138241` in mesh
millimetres (about 0.019 and 0.138 micrometres).  Thus the active build preserves
the valid multigrid hierarchy but does not claim exact radius on newly inserted
axial rings.  Achieving both requires repairing the coarse axial sweep before
TetGen—most likely station redistribution around high curvature without adding
near-coincident PLC rings—not changing the refinement projector alone.

Reproduction script:

`/Users/alhilaba/Documents/ promesh/work/test_ring_projection_audit.lua`

Validated inspection meshes from the final safe run are currently named:

- `/Users/alhilaba/Documents/ promesh/work/ring_radius_fixed_refined_1.ugx`
- `/Users/alhilaba/Documents/ promesh/work/ring_radius_fixed_refined_2.ugx`

The name reflects the experiment series; these two files contain the final
hierarchy-safe rule described above.

## W-shaped axial rings: artifact trace and final repair

The 19 September reference artifacts were inspected directly rather than
inferring their topology from the ProMesh display:

- `4th_trace_third.ugx` was born/modified at 2026-09-19 15:29:29;
- `4th_trace_third_refined_1.ugx` was born at 15:29:36 and modified at 15:29:37;
- the coarse grid contains 91,233 Inter tetrahedra, 9,248 Membrane
  hexahedra, and a Lumen made from 4,624 hexahedra plus 4,624 prisms;
- it contains no pyramid volumes;
- the refined grid preserves those region/type proportions exactly (eight
  children per parent volume).

All 4,624 old Apical axial edges and all 9,248 old Basolateral axial edges
have their child vertex at the exact Cartesian parent midpoint.  Thus the old
artifact also used axial linear interpolation; its geometry merely made the
inward sagitta too small to notice.  Its transverse edges were projected.

On the current trace, full surface projection removed the W shape and reduced
the Apical/Basolateral radius residual to numerical noise, but exposed 18
invalid cells in two tight intervals (`t=0.82132..0.822135` and
`t=0.823766..0.824583`).  Those intervals coincide with adjacent 0.5-um SWC
samples turning by 22.7 and 25.5 degrees.  The old midpoint rule hid this
coarse-centerline defect by pulling child rings inward.

The verified repair is therefore:

1. use exact projector positions for axial as well as transverse children;
2. retain the deterministic rotation-minimizing frame (1,024 fixed transport
   steps); and
3. soften input control points above 10 degrees for 50 passes before cubic
   interpolation, instead of the former 25-degree/35-pass defaults.

Validation with anisotropy 1, O-grid 10, BoxRefs 4 and two global conforming
refinements produced:

- coarse: 29,211 vertices and 73,591 volumes;
- level 1: zero corrected or indeterminate/degenerate volumes;
- level 2: zero corrected or indeterminate/degenerate volumes;
- level-1 maximum surface-radius residual: `4.15097e-10` mm;
- level-2 maximum surface-radius residual: `4.15097e-10` mm.

Inspection meshes are `work/ring_smooth10_l2.ugx`,
`work/ring_smooth10_l2_refined_1.ugx`, and
`work/ring_smooth10_l2_refined_2.ugx`.
