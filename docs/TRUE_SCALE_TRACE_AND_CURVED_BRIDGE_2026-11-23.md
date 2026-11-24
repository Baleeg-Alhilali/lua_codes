# True-scale traces and curved hex Inter bridge — 2026-11-23

Branch: `feature/true-scale-trace-repair` (Model). The neuro_collection plugin
directory is not a git repository; its changed source is checkpointed under
`codex_history/checkpoints/2026-11-23_curved_hex_bridge/`.

Status: **geometry pass** for one true-scale trace (fourth trace, first third)
through r2. No flow solve has been run on these meshes.

## 1. Trace scale error

The SWC traces are in image pixels; the user confirmed 1 x 1 x 1 um voxels.
The pipeline default `-centerlineScale 24` treated every pixel as 24 um.
`3rd_trace_243_third.swc` (349 um raw, 290 um after the corner smoother) was
therefore meshed as a 6.96 mm tube (7.52 mm with OuterWall terminal
extensions). The full `3rd_trace_243.swc` is 1318.8 um. All earlier flow
results on `third_trace_third_hex_sim_ready.ugx` are 24x too long; Poiseuille
pressure drop and Apical area scale with that length.

At 1 um voxels a 10 um lumen + 10 um membrane tube does not fit the manual
traces: the third-trace segment approaches itself to 23.1 um (needs 44.4 um for
the preflight, more with the Inter shell) and the clicked corners have 2-4 um
circumradii.

## 2. New tools and options

- `mesh_generation/tools/repair_trace_clearance.py`: moves the centreline only
  where it violates non-local clearance or minimum bend radius; endpoints fixed,
  optional frozen wall-normal tails (`--frozen-head/--frozen-tail`), carries
  `# UG4_` OuterWall box metadata, reports length change and displacement.
  SHA-256 `0d81b018364a1926192ceec28fa6f50b966d9f4ed22894579a67508138aa002e`.
- Pipeline `-skipSmoothing 1` keeps the supplied trace; a warning is printed
  whenever `-centerlineScale` is not 1.
  SHA-256 `210cc92fdedbf07fa78a1a1a2e87af4d4ff6a1f2c4df60c008e64ca0ca4495d5`.

Workflow: repair coil -> pipeline run (writes `_outerwall_extended.swc`) ->
repair the extended trace with 40 um frozen tails -> final run with
`-skipSmoothing 1 -centerlineScale 1 -extendTerminals 0`.

## 3. Generator defect and fix (neurites_from_swc.cpp)

Projected refinement inverted the hex-Inter pyramid/tetrahedron bridge at
true-scale curvature (and caused the Nov 22 ogrid 8/10 sweep crashes). The
coarse bridge pyramids were 0.61-0.65 um tall with ~11 um bases on a 21.5 um
ring; the chord sagitta of the curved base is ~0.72 um. The bridge apex was a
linear centroid whose parameters copied corner `v[0]`'s angle and axial value,
and the bridge volumes were linear while their base face was projected, so the
refined base centre moved past the apex.

Fix: the apex is placed on the shell map at the cell-averaged axial position,
circular-mean angle and mean radius, and the whole bridge (volumes, faces,
edges, vertices) uses projSH subset 0. Only the far TetGen tetrahedra remain
linear. A preliminary edit that only made bridge edges/faces linear had no
effect (17/2/3 degenerate cells unchanged) and was reverted.

Hashes: pre-edit `3e1f2056f7650889e30986b850373f8f9d9c2f1578549ff2c07b7ea65207094f`
(copy in `rollback_archive/2026-11-23_curved_hex_bridge_pre_edit/`, with the
pre-edit plugin dylib); post-edit
`ae2b1a2a3c914e813c71e4abb1cc4576502c63220090240f407f0958c28e85b9`.

## 4. Validation (ogrid 12, BoxRefs 3, hex Inter, shell 3 um, r10/t10 um)

| Case | Min bend / clearance (um) | Before fix | After fix |
|---|---|---|---|
| 4th third, repair defaults | 22.8 / 52.9 | 17 degenerate at r1 | 1 at r1 |
| 4th third, bend 26 | 26.1 / 55.2 | 2 at r2 | 4 at r2 |
| **4th third, bend 28, clearance 55** | 28.0 / 56.8 | 3 at r1 | **r0/r1/r2 PASS** |
| 3rd third, variant A (bend 26) | 26.1 / 52.6 | r2 fail | 6 at r2 |
| 3rd third, variant B (clearance 60) | 23.0 / 60.5 | r1 fail | 2 at r2 |

Passing mesh: `runs/true_scale_trace/fourth/rc28c55_bridge.ugx` (r0 10,266
vertices / 31,346 volumes); minimum hex scaled Jacobian 0.934 / 0.975 / 0.988
at r0/r1/r2; zero inverted or degenerate cells. Its centreline moved at most
8.5 um (mean 2.0 um) from the manual trace; length 343 um vs 372 um resampled.

12 of 13 remaining failures lie on the inner side of bends with local radius
25-30 um: with a 23 um shell outer radius the inner side is compressed to
about 10-20 % of its axial length and the cylindrical map approaches a
singularity. Practical requirement: minimum bend radius about 28-30 um.

24x-scale regression (third trace, numRefs 1): o6_b4_a4_hit10 and
o10_b3_a4_hit20 still pass; previously crashing o8_b2_a4_hit10 and
o10_b2_a4_hit10 now pass.

## 5. Known limits

- The third-trace segment cannot reach a 28 um bend radius with the current
  repair: without length preservation the coil contracts (163 um at bend 28).
  A length-preserving repair is required.
- The bridge fix is validated serially in the generator only; the solver's
  in-memory MPI refinement of these meshes, flow solves, and lumen radial
  resolution (`-lumenRings`) remain to be done.

## 6. Later the same day: corrections and robustness work in progress

**Correction to section 4.** `validate_ugx_mesh.py` checks tetrahedra,
hexahedra and prisms but **not pyramids**, and it has no minimum-quality
threshold for slivers. A new in-memory checker
(`check_local_nephron_refinement`, below) found that the "passing" fourth-trace
mesh `rc28c55_bridge` has, after one refinement, a bridge pyramid with corner
scaled Jacobian -8.9e-5 and a tetrahedron with 2e-5. Its r0-r2 pass is
therefore a topology/orientation pass only, not a quality pass.

Changes made (all on this workstation, backups and hashes in
`rollback_archive/2026-11-23_radial_fade_pre_edit/`):

1. **Radial curved-to-linear fade (UG4 core, NeuriteProjector).** Optional
   `set_radial_fade(start, end)`; a refinement node with averaged radial
   parameter r in (start, end) is placed at x_lin + w (x_curved - x_lin), w a
   C1 smoothstep. Off by default; serialized as class version 1 (old UGX files
   load with the fade off — verified on `rc28c55_bridge.ugx`).
2. **Generator fade band.** `configure_nephron_inter_fade` / pipeline
   `-interFadeDistance` (um, default 10): Inter tetrahedra within the band get
   shell-map parameters and use the projector; tetrahedra spanning two tube
   stretches or wrapping the axis stay linear.
3. **Canonical TetGen input order** (boundary faces sorted by centroid, points
   by coordinates, canonical face rotation).
4. **Opt-in double-precision TetGen input** (UG4 core wrapper,
   `UG_TETGEN_DOUBLE_PRECISION`); the nephron importer enables it unless the
   variable is set. Fourth-trace coarse build: float 8/12 successes, double
   12/12. TetGen output still varies run to run with identical input
   (divergence inside TetGen segment recovery).
5. **`check_local_nephron_refinement`** (plugin, Lua-registered): hanging-node
   refinement near the nephron (projSH subset 0 + buffer), in a sphere, or
   globally; reports leaf counts, per-type minimum corner scaled Jacobian
   (including pyramids), inverted cells and hanging-vertex offsets.

Open issues found by the checker (fourth trace, ogrid 12):

| Setup | Worst pyramid J | Worst tetrahedron J | Generator r2 |
|---|---|---|---|
| 3 um shell, 4 layers, fade 10 | -8.9e-5 | 5e-8 (TetGen sliver on bridge surface) | 18 degenerate |
| 6 um shell, 2 layers, fade 10 | 0.013 | 6e-5 | 1 degenerate |
| 6 um shell, 2 layers, fade 20 | 0.013 | 5e-5 | 1 degenerate |

- Bridge cells are flat (thickness << 12 um circumferential size); a thicker
  bridge helps but does not remove the problem. Proposed fix: all-hex shell
  plus an outward pyramid cap on each outer face.
- TetGen leaves slivers whose four vertices lie on the curved bridge surface;
  the wrapper requests only a radius-edge bound (`pqq<q>Y`), no dihedral bound.
- The fade-band inverse map (station interpolation) is accurate only to
  0.5-0.9 um mean, 1.2-6.8 um max; it needs a Newton correction against the
  projector.

None of these meshes is accepted for simulation yet.

## 7. Outward-cap shell redesign and final r0-r2 checks (late 2026-11-23)

Scope rule from the user: the mesh must support global refinement well beyond
r2 (r5+ on a larger machine), but nothing beyond r2 is run on this workstation.
The evidence for deeper levels is therefore the *trend* r0 -> r1 -> r2.

### Changes (neuro_collection generator unless noted)

1. **All-hex shell with outward pyramid caps.** Every structured near-Inter
   layer is a hexahedron (uniform layers). Each outermost quad carries an
   outward pyramid whose apex lies on the tube map at 0.5 x the shortest base
   edge; the four cap triangles form the TetGen boundary. Replaces the
   centre-apex pyramid/tetrahedron bridge (flat cells) and its
   address-dependent diagonal choice. Prism Inter path unchanged.
2. **Interior shell/cap vertices placed at the projector image of their
   parameters** before TetGen (terminal-plane vertices keep their snapped wall
   positions). Their geometric construction differed from the map by
   0.11 um mean / 2.5 um max; refinement nodes come from the map, so that
   mismatch compounded into slivers. This was present in the hex shell before
   this session.
3. Fade band: safeguarded Newton correction of the inverse map (mean residual
   ~1-4 nm), demotion of inconsistent or poorly conditioned vertices, closure
   rule (a projected tetrahedron projects all its edges and faces).
4. Opt-in TetGen minimum dihedral (`UG_TETGEN_MIN_DIHEDRAL`, UG4 core
   wrapper). Not used: 10-20 degrees did not raise the boundary-pinned coarse
   sliver quality and 20 degrees quadrupled the cell count.
5. **Fade default changed to off** (`-interFadeDistance 0`). With caps, the
   fourth trace keeps tetrahedron quality 0.0195 at r1 and r2 without fade;
   with a 10 um fade it settled at 1e-4 (20 um: 5e-6 at r2). Tube-coordinate
   interpolation offsets (~0.4 um on 10 um coarse edges) exceed the thickness
   of the coarse TetGen slivers adjacent to the caps. The fade code remains
   available for meshes with well-shaped tetrahedra in the band.
6. Checker mode 3 (linear refinement, projector disabled) and level-0 /
   far-field statistics. Linear refinement of the same coarse mesh keeps every
   cell type's quality exactly (tet 0.038 at r0, r1, r2), which localized the
   earlier degradation to projected placement.

Source checkpoint: `codex_history/checkpoints/2026-11-23_caps_fade/`
(SHA256SUMS inside). Checker driver: `mesh_generation/check_refinement_quality.lua`.

### Final r0-r2 results (ogrid 12, BoxRefs 3, hex Inter, fade off)

Minimum corner scaled Jacobian (generator validation + in-memory checker to r2):

| Case | Min bend (um) | Result | tet r1/r2 | hex r1/r2 | pyramid r1/r2 |
|---|---:|---|---|---|---|
| **4th third, bend 28, shell 3 um** | 28.0 | **r0-r2 PASS** | 0.0195 / 0.0195 | 0.949 / 0.964 | 0.129 / 0.148 |
| 4th third, shell 6 um, 2 layers | 28.7 | r1 inverted cap | - | 0.43 (r0 0.39) | -0.005 |
| 4th third, repair defaults | 22.8 | generator segfault | - | - | - |
| 3rd third, variant A | 26.1 | r1 degenerate (4) | - | - | - |
| 3rd third, variant B | 23.0 | r1 degenerate (3) | - | - | - |
| 24x regressions (4 configs) | n/a | r0-r1 PASS | o6: 0.0003 / 0.0003 | 0.73 / 0.76 | 0.016 / 0.071 |

The passing true-scale case keeps (or improves) every minimum quality from r1
to r2, which is the property needed for deeper refinement; r3+ is untested.

Failures are geometric: the structured region (outer radius 20 um + shell +
cap apex at ~29 um for the 3 um shell) approaches or exceeds the local bend
radius on the inner side of bends, where the tube map folds. A thicker shell
makes this worse. Needed next: curvature-limited cap height (like the
existing curvature-limited shell thickness), and a repair-tool default that
derives the minimum bend radius from outer radius + shell + cap.

## 8. Curvature-limited caps (rejected) and default workflow

- Tried limiting cap height so the apex stays within 0.75 x local bend radius
  (floor 0.15 x shortest base edge). Rejected and reverted: lowered caps are
  flat pyramids again — the previously passing fourth-trace case inverted at
  r1, and the third-trace cases only improved from 4/3 to 2/2 degenerate cells.
  The generator is back at the checkpointed version
  (`c61a0c304dc4904824f000e28436eb7581d592cc770aa3ff2f2ce446c9e419db`).
- `repair_trace_clearance.py` default minimum bend radius is now
  `(outer + inter shell)/0.82` (28.0 um for 20 + 3 um), calibrated on the r0-r2
  results above.

End-to-end default workflow (repair -> extend -> repair with 40 um frozen
tails -> generate r0-r2 -> in-memory check to r2), `runs/true_scale_trace/default_workflow/`:

| Trace | Length raw -> repaired | Max / mean move | Generation | tet r1/r2 | hex r1/r2 | pyramid r1/r2 |
|---|---|---|---|---|---|---|
| 4th third | 380.5 -> 316.6 um | 21.7 / 4.0 um | r0-r2 PASS | 0.0089 / 0.0052 | 0.949 / 0.963 | 0.142 / 0.147 |
| 3rd third | 349.1 -> 149.2 um | 53.8 / 22.3 um | r0-r2 PASS | 5.4e-4 / 1.7e-5 | 0.857 / 0.747 | 0.006 / 0.028 |

The fourth trace is a usable true-scale mesh (stable hex/prism/pyramid
quality; tetrahedra decline mildly). The third-trace segment cannot hold a
20 um tube at a 28 um bend radius without the coil contracting to less than
half its length, and its refined quality still declines; it needs either a
length-preserving repair or a check of the tracing/radius against the image.
The earlier hand-tuned fourth-trace repair (bend 28 / clearance 55, max move
8.5 um, 343 um) remains the most faithful passing case.
