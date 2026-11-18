# Scale-adaptive nephron geometry

Date: 2026-10-08

## Goal

Allow the same mesh pipeline to accept different lumen radii and membrane
thicknesses, including membrane thicknesses of 10--15 um, without relying on
absolute tolerances that only work for the original 1.0 um lumen radius and
0.2 um membrane.

The production entry points covered by this change are:

- `work/outerwall_portal/single.lua`
- `work/outerwall_portal/multiple.lua`
- `work/outerwall_portal/neurites_from_swc.cpp`
- `extend_swc_terminals_to_outerwall.py`
- `quadratic_interpolate_swc.py`

No UG4 core or NavierStokes source was changed.

## Root cause found

The failures had two separate causes.

1. Some mesh-cleanup and geometric comparison tolerances in the Neuro
   Collection changes were absolute mesh-coordinate values. They were valid
   for the reference radius but did not retain the same meaning after a large
   change in tube or box size.
2. A 10 um lumen radius plus a 10--15 um membrane does not physically fit the
   unscaled third-trace centerline. The requested outer radius is 20--25 um,
   while the measured minimum discrete centerline curvature radius is only
   5.115485 um. The swept tube self-intersects before TetGen is called. This is
   not repairable by relaxing a numerical tolerance.

The old behavior allowed the invalid surface to reach TetGen, where it could
segfault. The pipeline now measures this condition and stops with a calculated
minimum `-centerlineScale`.

## Relative controls now used

| Operation | Previous behavior | Scale-adaptive behavior |
|---|---|---|
| terminal extension sampling | fixed 0.5 um | `max(resampleSpacing * centerlineScale, 0.25 * outerRadius)` |
| non-terminal box padding | fixed 10 um | `max(10 * centerlineScale, 2 * outerRadius)` |
| terminal box padding | already relative, 25 radii | remains `terminalPaddingFactor * largestOuterRadius` |
| terminal snap numerical floor | implicit/fixed floor | machine epsilon times OuterWall diagonal |
| terminal snap acceptance | geometry-dependent constant | 2% of OuterWall diagonal |
| coincident projector vertices | fixed decimal quantization | quantization and comparison based on grid diagonal |
| OuterWall plane detection | fixed positional floor | machine epsilon and `1e-10 * boxDiagonal` |
| recorded-volume containment | fixed `1e-9` mesh units | `1e-6` of the recorded element scale |
| duplicate mesh vertices | fixed physical tolerance | `1e-4 * minimumOuterRadius` |
| projector-radius audit | fixed absolute residual | `1e-8 * requestedRadius` |

The remaining small constants in axial/radial projector coordinates are
dimensionless parameter-space comparisons. Angle limits, interpolation
fractions, quality ratios, and refinement marks are also dimensionless and
must not be multiplied by the physical radius.

## New geometry parameters

- `-radius R`: lumen radius in um.
- `-membraneThickness T`: physical membrane thickness in um.
- `-centerlineScale S`: multiplies SWC x/y/z while leaving the requested tube
  radii unchanged. Use this to enlarge the complete nephron layout when the
  larger cross-section cannot fit its original bends or clearances.
- `-terminalSpacing H`: optional terminal-extension sampling distance in um.
- `-padding P`: optional OuterWall clearance in um.
- `-paddingX PX`, `-paddingY PY`, `-paddingZ PZ`: independent symmetric
  clearance from the nephron outer surface to the two box faces normal to each
  axis. Each defaults to `-padding`.
- `-skipGeometryPreflight 1`: diagnostic escape hatch only; it should not be
  used for production meshes.

Intermediate ready and extended SWCs are now named from `-out`. Parameter
sweeps therefore no longer overwrite the shared source-side `_mesh_ready.swc`.

The terminal curvature floor is applied after the axis values are read. With
the default `-terminalPaddingFactor 25`, each clearance is at least 25 outer
radii. To request a smaller exact axis clearance, reduce that factor as well;
use `0` only after checking the terminal transition curvature and final mesh.

Example with independent box clearance:

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex work/outerwall_portal/single.lua -- \
  -swc 'Nephron traces/3rd_trace_243_third.swc' \
  -radius 1 -membraneThickness 0.2 \
  -paddingX 40 -paddingY 50 -paddingZ 70 \
  -out runs/axis_padding_test.ugx
```

This produces 40 um clearance on each X side, 50 um on each Y side, and 70 um
on each Z side.

## Preflight rules

For outer radius `R_outer = radius + membraneThickness`, the default check
requires:

```text
R_outer < 0.9 * minimum centerline curvature radius
2 * R_outer < 0.9 * minimum non-local centerline clearance
```

The 0.9 factor is a dimensionless safety margin. It leaves room for polygonal
cross-sections and interpolation error instead of treating exact tangency as a
valid mesh.

## Verified tests

### Reference-size regression

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex work/outerwall_portal/single.lua -- \
  -swc 'Nephron traces/3rd_trace_243_third.swc' \
  -radius 1 -membraneThickness 0.2 -numRefs 1 \
  -out runs/scale_regression_r1_t0p2.ugx
```

Result: PASS.

- coarse vertices: 29,211
- coarse volumes: 73,591
- refined indeterminate/degenerate volumes: 0
- maximum Apical/Basolateral radius residual: `4.15097e-10`

### Large variable geometry

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex work/outerwall_portal/single.lua -- \
  -swc 'Nephron traces/3rd_trace_243_third.swc' \
  -radius 10 -membraneThickness 15 -centerlineScale 6 -numRefs 1 \
  -out runs/scale_valid_r10_t15_s6_final.ugx
```

Result: PASS.

- requested outer radius: 25 um
- measured scaled minimum curvature radius: 30.692911 um
- measured non-local centerline clearance: 86.282704 um
- coarse vertices: 13,187
- coarse volumes: 33,410
- refined indeterminate/degenerate volumes: 0
- maximum Apical/Basolateral radius residual: `1.10018e-9`

The corresponding unscaled request is rejected before meshing and reports:

```text
Requested outer radius 25 um is incompatible with centerline curvature
(minimum radius 5.11549 um). Use -centerlineScale > 5.43014.
```

## Source fingerprints

```text
239f29de6a4a400264ae94655022b681d7ad2e041c2b12364c2a7548901e446d  work/outerwall_portal/neurites_from_swc.cpp
a782831a21cc0161e1f59d47cd811018aa6d11fc4828a33cc5441cfc3711b05a  work/outerwall_portal/single.lua
53168d97c2db1d871285ae1a8c58c1200bddd83fa0870a28cd30143fd39cc7c2  work/outerwall_portal/multiple.lua
579be6b849091b234b20ec47bdcfd7e15d50091cf605a47bdda02cb9211c1088  runs/scale_valid_r10_t15_s6_final.ugx
```
