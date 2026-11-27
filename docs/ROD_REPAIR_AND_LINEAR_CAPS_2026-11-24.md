# Length-preserving trace repair and linear outward caps — 2026-11-24

Branch `claude/mesh-generation-d62595` (Model, from `aff59b4`). UG4 tree
`/Users/alhilaba/UG4_promesh_ogrid` (ugcore `46b9b57a8` + overlay). Serial runs on
the 36 GB workstation, nothing above refinement level 2.

Status: **geometry pass r0-r2 with stable r1 -> r2 quality** for the third
trace (segment and full trace). No flow or transport solve has been run on
these meshes. r3+ is not run here; the r1 -> r2 trend is the evidence.

Follows `TRUE_SCALE_TRACE_AND_CURVED_BRIDGE_2026-11-23.md`.

## 1. Length-preserving repair (`mesh_generation/tools/repair_trace_rod.py`)

The earlier relaxation (`repair_trace_clearance.py`) leaves the length free,
and curvature relaxation acts as curve-shortening flow: the third-trace
segment contracted 349 -> 149 um, the full fourth trace 1134 -> 942 um.

The new tool models the centreline as an **inextensible rod**: N segments of
fixed length h (2.5 um), minimized with L-BFGS (PyTorch, double):

- inextensibility `((|x_{i+1}-x_i| - h)/h)^2` (weight 1e4),
- bend: three-point circumcurvature at h and 2h above 1/R* (squared excess),
- clearance: every pair farther apart than 80 um of arc closer than C*,
- weak fidelity to the 2 um Gaussian-denoised trace at the same arc length
  (arc length is material for an inextensible rod), weak fairing, end anchor,
- a homotopy ramps R*, C* from the input values to the targets (+3 % margin);
  the constraint weight is then raised x4 per stage until the 5 um output
  passes the same audit as the Lua preflight.

The excess length has to go somewhere, so coils expand outward and lift out
of their plane (see the PNGs). Target length = denoised length (raw minus
voxel-click noise, ~3 %). Endpoints fixed by default (`--free-ends` anchors
them softly), `--frozen-head/--frozen-tail` keep the OuterWall tails rigid and
carry the `# UG4_` box metadata. Several inputs are repaired **jointly**
(clearance also enforced between traces, without arc exclusion).

Driver: `mesh_generation/tools/rod_repair_workflow.sh <raw.swc> <out_dir> <prefix>`
(coil repair -> pipeline pass writing `_outerwall_extended.swc`, padding 100 um
-> frozen-tail repair -> mesh r0-r2 -> optional `CHECK=1` in-memory check).

### Bend target 32 um, not 28 um

With the 28 um limit ((outer + shell)/0.82) the rod repair of
`3rd_trace_243_third` (349 -> 339 um) generated r0, r1 and failed at r2 with one
degenerate cap pyramid on the inner side of the bend at the end of a 40 um arc
held exactly at R = 28.9 um. The free-length repair only touched the limit at
isolated points; the rod holds long arcs on it. Default is now
(outer + shell)/0.72 = **32 um**. 35 um was infeasible with fixed ends.

### Results (1 um voxels, clearance 55 um, bend 32 um)

| Trace | Raw -> repaired | Kept | Min bend / clearance | Max / mean move |
|---|---|---:|---|---|
| 3rd_trace_243_third | 349.1 -> 343.5 um | 98.4 % | 32.4 / 55.4 | 24.8 / 11.0 um |
| 3rd_trace_243 (full) | 1318.8 -> 1293.5 um | 98.1 % | 32.1 / 55.4 | 25.0 / 8.0 um |
| 4th_trace_167 (full) | 1134.1 -> 1114.7 um | 98.3 % | 32.4 / 56.3 | 23.4 / 7.4 um |
| 4th_trace_third | 380.5 -> 370.2 um | 97.3 % | 32.2 / 61.4 | 9.1 / 3.3 um |

Previous free-length repair: 149 um (3rd third, 43 %), 942 um (4th full, 83 %).
Each run takes 15-70 s. Files: `Nephron traces/true_scale_repaired/rod_2026-11-24/`.

### Joint repair of the third and fourth traces (for one shared Inter mesh)

The raw full traces come within 18.3 um of each other; independent repairs
leave 4.9 um (about 250 points closer than 55 um). Both traces end at the
image edge (y ~ 110 um) with terminals 25.5 um apart, so fixed endpoints make
the joint problem infeasible. `--free-ends` joint repair:

| Trace | Raw -> repaired | Max / mean move |
|---|---|---|
| 3rd_trace_243 | 1318.8 -> 1286.8 um (97.6 %) | 29.6 / 11.4 um |
| 4th_trace_167 | 1134.1 -> 1114.6 um (98.3 %) | 28.8 / 10.3 um |

Min bend 32.8 um, min self-clearance 56.1 um, **min distance between the
traces 55.1 um**. Files `*_joint_with_*_coil.swc`, `joint_3rd_4th_coil.json`.
Not yet meshed: the shared mesh needs the multi-nephron pipeline, terminal
extension of both traces to the OuterWall, and a joint frozen-tail repair.

## 2. r1 -> r2 tetrahedron decline: cause and fix

Full fourth trace (`mesh_generation/meshes/2026-11-24_fourth_pad100/fourth_trace_full_pad100.ugx`),
checker mode 2: tet 0.0140 (r1) -> 0.00126 (r2); hex/prism/pyramid stable.

### Diagnosis

The checker (`check_local_nephron_refinement`) now prints, at every level, the
minimum J of tetrahedra with and without a projected (projSH 0) edge, and the
worst tetrahedron's ancestry (vertices with parameters, edge/face projSH).

- r1 worst: a linear TetGen tetrahedron in the gap between two passes of the
  coil (corner axial parameters 0.37 and 0.50, ~160 um apart along the path,
  10 um apart in space), spanning cap vertices of both passes.
- r2 worst: a linear tetrahedron in the valley between two circumferentially
  adjacent caps, with two apexes, their shared base corner and a far vertex.
  Ancestry J: 0.076 (r0) -> 0.027 (r1) -> 0.0013 (r2).

In both cases the projected edges are the **slanted cap edges** (base corner
-> apex). The projector interpolates (axial, angle, radius) along them, so in
the cross-section the edge is a spiral: its midpoint lies 0.33 um (3.8 % of
the edge) off the chord at r1 and 2 % at r2 (measured). The neighbouring
tetrahedra are linear and have small corner angles at cap vertices. The
corner child keeps the parent shape, so each level's offset removes about
the same fraction of thickness: losses 0.049, 0.026, ... sum to ~0.1 > 0.076,
i.e. that tetrahedron would invert near r3. Linear refinement (mode 3) keeps
all qualities exactly, which confirms projected placement as the cause.

### Fix (neurites_from_swc.cpp)

The cap is only the interface to TetGen, so its curvature is not needed:
cap pyramids, their triangles and slanted edges are now projSH 1 (linear);
the base quad (shared with the hex shell) and its edges stay on the tube map;
the apex keeps its map position. `UG_NEPHRON_PROJECTED_CAPS=1` restores the
fully projected caps. The checker's mode 4 applies the same reassignment to an
existing UGX (diagnostic, identical coarse mesh, no TetGen variation).

Same coarse mesh, full fourth trace:

| Mode | tet r1 / r2 | hex r1 / r2 | pyramid r1 / r2 |
|---|---|---|---|
| 2, projected caps (old) | 0.0140 / 0.00126 | 0.954 / 0.952 | 0.134 / 0.147 |
| 4, linear caps | **0.0213 / 0.0211** | 0.954 / 0.952 | 0.081 / 0.073 |

Tetrahedra now refine almost exactly (the minimum is the coarse TetGen value).
Linear caps lower the pyramid minimum (the curved base bulges toward a linear
apex), but this loss shrinks with each level (decrements 0.008 after r1) toward a
positive limit because the cap has real height (0.5 x shortest base edge),
unlike the earlier flat bridge cells.

## 3. Third trace meshes with linear caps (generator r0-r2 + checker mode 2)

Settings: radius 10, membrane 10, ogrid 12, BoxRefs 3, hex Inter, shell 3 um,
`-skipSmoothing 1 -centerlineScale 1 -extendTerminals 0`, padding 100 um.

| Mesh | Coarse V / volumes | Generation | tet r1 / r2 | hex r1 / r2 | prism r1 / r2 | pyramid r1 / r2 |
|---|---|---|---|---|---|---|
| 3rd third (343 um coil, 627 um with tails) | 9,652 / 27,888 | r0-r2 PASS | 0.0179 / 0.0179 | 0.967 / 0.968 | 0.469 / 0.471 | 0.157 / 0.139 |
| 3rd full (1294 um coil, 1585 um with tails) | 23,982 / 73,379 | r0-r2 PASS | 0.0266 / 0.0266 | 0.966 / 0.978 | 0.469 / 0.470 | 0.146 / 0.136 |
| 3rd third, projected caps (before fix) | 9,652 / 27,888 | r0-r2 PASS | 0.0077 / 0.0039 | 0.967 / 0.968 | 0.469 / 0.471 | 0.197 / 0.147 |

Zero inverted or degenerate cells at every level; the worst tetrahedra of the
full trace are now pure far-field TetGen cells (radial = 0) that refine without
loss. r2 leaf counts: 1.80 M (third) and 4.75 M volumes (full).

Outputs (not in git): `runs/rod_repair_2026-11-24/lc2/` (third_b32_lc,
third_full/), earlier variants in `runs/rod_repair_2026-11-24/`. The folder
`runs/rod_repair_2026-11-24/lc/` is an invalid run (zsh passed the mesh
flags as one word; trace scaled x24) and can be deleted.
Coarse full-third mesh for solvers: `mesh_generation/meshes/2026-11-25_third_pad100/`.

## 4. Known limits / next steps

- The full fourth trace has not been regenerated with linear caps (checked
  only via mode 4 on the existing coarse mesh); TetGen output varies run to run.
- Pyramids decline mildly r1 -> r2 with linear caps (about 7-11 %); r3+ untested.
- A throwaway extension pass (`-numRefs 0`) segfaulted inside TetGen once on
  the rod-repaired coil; the extended SWC is written before TetGen, and the
  final runs did not crash. TetGen robustness is not resolved.
- The joint third+fourth mesh is not built yet (see section 1).
- The 32 um bend radius was calibrated on these traces only.

## 5. Changed files and hashes

- `mesh_generation/tools/repair_trace_rod.py`, `mesh_generation/tools/rod_repair_workflow.sh` (new).
- `plugins/neuro_collection/grid_generation/neurites_from_swc.cpp`:
  `c61a0c30...` -> `55402eebee844c48417a62f237b8db3c3a2d4c6f7ff5aac8e15f172aff5b3fab`
  (linear caps, checker lineage + mode 4). Pre-edit copies with SHA256SUMS:
  `rollback_archive/2026-11-24_mesh_tet_decline_pre_edit/` (main checkout).
  Overlay updated (`ug4_build_overlay/`); all neuro_collection hashes verified
  against the installed tree. Rebuilt with `make neuro_collection` only.

## 6. Two traces in one shared Inter domain (later the same day)

`mesh_generation/swc_to_neurite_projector_multiple_nephrons.lua` was updated
to the true-scale pipeline:

- options from the single pipeline: `-skipSmoothing`, `-interMesh hex|prism`,
  `-prismInterThickness`, `-interFadeDistance`, the geometry-policy options,
  structural validation of the coarse and every refined file (`-validate`);
  defaults aligned (radius 10 / membrane 10 um, `bendChordErrorFactor` 0.01,
  `resampleSpacing` 5, `terminalSpacing` 20, `terminalPaddingFactor` 0 — the
  old 25 x outer radius gave a 500 um minimum box padding);
- `-mergedSwc <file>`: a pre-merged multi-root SWC in um (the joint rod repair
  of this script's own `_outerwall_extended.swc`); no smoothing, merging or
  extension, `# UG4_` box metadata reused;
- a preflight of non-local clearance within each trace and the distance
  between traces (2 x outer / 0.9);
- `-minPortSeparation` (default 3 x outer = 60 um), passed to
  `extend_swc_terminals_to_outerwall.py --min-port-separation`: independent
  nearest-cell snapping put both outlets (both on ymin) in adjacent cells
  48.5 um apart; a later port on the same wall now moves to the nearest
  cell centre at least that far from earlier ports.

`repair_trace_rod.py` reads a merged multi-root SWC as several curves and
writes it back as one multi-root file (joint frozen-tail repair).

**C++ (`neurites_from_swc.cpp`).** The structured near-Inter shell builder
(`add_straight_prismatic_inter_to_version2_nephron`) grouped Basolateral rings
by axial parameter only; with two nephrons both run over [0, 1], so the rings
mixed ("Could not propagate one-to-one Basolateral ring ordering"). Stations
are now keyed by (nephron, axial); ring-order propagation, terminal ramps,
curvature limits, slope limits, the sweep and the Basolateral face assignment
stay within one tube. The TetGen fill, portal rebuild and caps were already
global. The fade band throws for more than one nephron. The multi-nephron entry
also enables double-precision TetGen input like the single one. Single-trace
regression: the third segment rebuilds identically (9,652 vertices / 27,888
volumes, same shell statistics). Hash
`300129f941927a4ff1e031556468c8e15ecbb6af0e0f0979cdf61e8709a64ae4`; backup of
the previous file and the multi script in
`rollback_archive/2026-11-24_multi_nephron_pre_edit/`.

**Workflow** `mesh_generation/tools/rod_repair_joint_workflow.sh`: joint coil
repair with `--free-ends` (the raw traces end 25.5 um apart at the image edge,
which makes fixed ends infeasible) -> multi script pass writing the merged
extended SWC -> joint frozen-tail repair -> `-mergedSwc` mesh -> checker.

**Bend 32 um failed for the joint case, 35 um passes.** With 32 um the fourth
trace (tree 2) had 4 inverted cells at r1 (linear cap pyramid -0.0069, tet
1e-4, hex 0.88) on the inner side of bends at R ~ 33 um; the same failure
appeared when that tree was meshed alone with the single pipeline, so it is
the trace geometry, not the multi-nephron code. With a 35 um target (now the
joint workflow default):

| Joint mesh (third + fourth) | value |
|---|---|
| traces | 1318.8 -> 1289.2 um (97.8 %), 1134.1 -> 1114.8 um (98.3 %), max move ~30 um |
| merged preflight | bend >= 35.3 um, self-clearance >= 55.9 um, between traces 56.2 um |
| coarse | 42,484 vertices, 122,420 volumes, validation PASS |
| checker r1 / r2 | tet 0.00222 / 0.00222, hex 0.971 / 0.979, prism 0.469 / 0.471, pyramid 0.191 / 0.139, 0 inverted |
| r2 size / memory | 7.93 M leaf volumes, 11.7 GB peak (checker) |

The worst tetrahedron is a linear far-field TetGen cell at the inlet wall that
refines without loss after r1. Pyramids lose ~27 % r1 -> r2 near a terminal
(axial < 0.01, where the shell thickness ramps); the third-trace meshes show
the same terminal-region trend (7-11 %). This is the main open quality item
before deeper refinement. Mesh: `mesh_generation/meshes/2026-11-25_third_fourth_joint/`.

Not done: generator r2 flat files for the joint mesh (r2 checked in memory),
flow/transport solves, and the 32 vs 35 um choice for the single-trace meshes
(the third-trace meshes used 32 um and passed).

## 7. Thicker Inter hex cells on the full fourth trace (2026-11-25)

The hex cell radial thickness is `-prismInterThickness / layers`, where the
layer count is clamp(ceil(`-interLayersPerSurfaceEdge` x thickness / surface
edge), `-interMinLayers`, `-interMaxLayers`); with the 3 um shell the minimum
of 4 layers gave 0.75 um cells. Lua and C++ both require >= 2 layers.
Scan (2 layers, same trace as `mesh_generation/meshes/2026-11-24_fourth_pad100`): 3 and 4 um
pass r0-r2 with 0 inverted cells; 4.5 um inverts 2 cap pyramids at r2; 5-7 um
segfault in TetGen. Meshes, commands and quality table:
`mesh_generation/meshes/2026-11-25_fourth_thick_hex/README.md`. The 2 x 2 um mesh is the thickest
valid one but its pyramid minimum drops 0.081 -> 0.020 r1 -> r2; the 2 x 1.5 um
mesh is stable. Geometry pass only; no solves.

Prism shell instead (`-interMesh prism -interRadialGrowth 0`, 2 layers, same
trace), 2026-11-25: coarse generation passes at 3, 4, 5, 6 and 8 um, but every
variant fails refinement with "Could not repair a negatively ordered refined
prism" (3, 4, 6 um at r1; 8 um at r2, after r1 prism J 0.061). The coarse
TetGen far field is also much worse than with hex: 1,715-3,825 tetrahedra
below J 0.05 (minimum 1e-6 to 2e-3) versus 23-34 with hex. The prism path has
no outward caps (2,520 triangular transition faces) and none of the hex
path's bend/portal thickness limits. Prisms do not fix the thin-cell problem;
the 2-layer hex meshes in `mesh_generation/meshes/2026-11-25_fourth_thick_hex/` remain the option.

## 8. Fourth trace with a 10 um, two-layer hex shell (2026-11-26)

A shell as thick as the membrane (10 um, 2 x 5 um hexes) needs bend radius
>= ~50 um (the builder caps the shell at 0.6 R - 20 um) and centreline
clearance >= ~70 um. The fourth trace was re-repaired at R 55 / clearance 75
(default homotopy reached R 51 only; `--ramp-steps 24 --lbfgs-iters 1000
--w-length 1e5` passed). Built with anisotropy 2 so that Inter/Membrane hex
aspect ratios stay below 10 (max 8.2 at r2). Validation PASS r0-r2, 0 inverted.
Details and tables: `mesh_generation/meshes/2026-11-26_fourth_shell10/README.md`.
`rod_repair_workflow.sh` now takes `INTER_OPTS`; new tool
`mesh_generation/tools/hex_aspect_ratio.py`. Geometry pass only.
Cap treatment follow-up: anisotropy 1.5 (squarer outer cap bases, 174
stations) raised the r2 pyramid minimum 0.022 -> 0.101 and the tet minimum
0.0054 -> 0.0221 with fewer volumes; finer TetGen near-field sizing made the
mesh 5-6x larger and the tets worse, and `-interQuality 1.4` had no effect.
Recommended mesh: `mesh_generation/meshes/2026-11-26_fourth_shell10/fourth_shell10_a15.ugx`.
