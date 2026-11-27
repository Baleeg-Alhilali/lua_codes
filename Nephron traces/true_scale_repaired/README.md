# True-scale repaired nephron traces (2026-11-23/24)

## Why the traces were changed

The SNT traces are in image pixels; the voxel size is 1 x 1 x 1 um. The
former pipeline default `-centerlineScale 24` treated every pixel as 24 um,
so meshes were 24 times too long. At true scale, a tube with a 10 um lumen,
10 um membrane (20 um outer radius) and 3 um structured Inter shell does not
fit the hand-clicked centrelines: they have corners with 2-4 um curvature
radius and parts of the trace pass within ~22 um of each other.

Raw audit (`mesh_generation/tools/audit_trace.py`, 1 um voxels, outer radius
20 um, shell 3 um; `audit_raw_*.json`):

| Trace | Length | Min bend | Points sharper than 28 um | Min self-distance |
|---|---:|---:|---:|---:|
| 4th_trace_167 (full) | 1134 um | 2.7 um | 82 | 21.9 um |
| 4th_trace_third | 381 um | 3.3 um | 36 | 47.8 um |
| 3rd_trace_243 (full) | 1319 um | 3.1 um | 120 | 22.1 um |
| 3rd_trace_243_third | 349 um | 3.8 um | 36 | 23.1 um |

Mesh limits (calibrated on r0-r2 meshes, `docs/TRUE_SCALE_TRACE_AND_CURVED_BRIDGE_2026-11-23.md`):
minimum bend radius ~(outer + shell)/0.82 = 28 um; non-local clearance
>= 2 (outer + shell)/0.9 = 51 um, 55 um used in practice.

## How they were shifted

`mesh_generation/tools/repair_trace_clearance.py` (position-based relaxation
on a 1 um resampled polyline; endpoints fixed):

1. push apart every pair of points farther than 80 um apart along the trace
   but closer than the clearance (displacements smoothed along the curve);
2. relax points whose circumradius (at 5 um sampling) is below the minimum
   bend radius toward their neighbours' midpoint;
3. weak global smoothing and uniform re-spacing; length is free and reported;
4. stop when both limits hold on the 5 um output sampling (same formulas as
   the Lua preflight).

The OuterWall terminal extensions are then added by the normal pipeline, and
the extended trace is repaired again with the straight wall-normal tails
frozen (`--frozen-head 40 --frozen-tail 40`), which keeps the OuterWall box
metadata valid; the mesh is built with `-skipSmoothing 1 -centerlineScale 1
-extendTerminals 0`.

Commands (from the Model root):

```sh
# 1. repair the coil (um in, um out)
python3 mesh_generation/tools/repair_trace_clearance.py "Nephron traces/4th_trace_167.swc" \
    --clearance 55 -o coil.swc --report coil.json
# 2. pipeline pass that writes <out>_outerwall_extended.swc (padding 100 um)
(cd mesh_generation && <UG4>/bin/ugshell -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- \
    -swc ../coil.swc -skipSmoothing 1 -centerlineScale 1 -padding 100 -radius 10 \
    -membraneThickness 10 -ogrid 12 -BoxRefs 3 -interMesh hex -prismInterThickness 3 \
    -numRefs 0 -validate 0 -out ../p1.ugx)
# 3. repair the extended trace with frozen tails (input in mm)
python3 mesh_generation/tools/repair_trace_clearance.py p1_outerwall_extended.swc --voxel-size 1000 \
    --frozen-head 40 --frozen-tail 40 --clearance 55 -o extended_repaired.swc --report ext.json
# 4. final mesh r0-r2 (see mesh_generation/meshes/2026-11-24_fourth_pad100/README.md)
# evaluate / plot
python3 mesh_generation/tools/audit_trace.py extended_repaired.swc --voxel-size 1
python3 mesh_generation/tools/plot_trace_repair.py "Nephron traces/4th_trace_167.swc" coil.swc -o cmp.png
```

The repair is path dependent: the default clearance (51 um) contracted the full
fourth trace to 577 um and failed meshing, while 55 um gave 942 um.

## Files and results

| File | Source trace | Settings | Length | Max / mean move* | Min bend / clearance | Mesh r0-r2 |
|---|---|---|---:|---:|---|---|
| `4th_trace_167_full_repaired_centreline.swc` | 4th_trace_167 | bend 28, clearance 55 | 942 um | 38.8 / 8.0 um | 28.3 / 56.6 um | — |
| `4th_trace_167_full_pad100_extended_repaired.swc` | above + tails | padding 100, frozen tails | 1206 um lumen path | 6.3 um (tails pass) | 28.1 / 56.6 um | PASS (used for flow) |
| `4th_trace_third_repaired_centreline.swc` | 4th_trace_third | bend 28, clearance 55 | 343 um | 8.5 / 2.0 um | 28.3 / 57.4 um | — |
| `4th_trace_third_pad200_extended_repaired.swc` | above + tails | padding 200 | — | — | 28.0 / 56.8 um | PASS |
| `3rd_trace_third_variantA_bend26_centreline.swc` | 3rd_trace_243_third | bend 26 | 302 um | 21.5 um | 26.1 / 52.6 um | fails (tight bend) |
| `3rd_trace_third_variantB_clear60_centreline.swc` | 3rd_trace_243_third | clearance 60 | 357 um | 30.9 um | 23.0 / 61.4 um | fails (tight bend) |

\* `*_report.json` measures moves on the 1 um working polyline; the PNG plots
measure the 5 um output points against the original and show smaller maxima
(e.g. 18.6 um for the full fourth trace).

The third-trace segment cannot reach a 28 um bend radius without the coil
contracting to less than half its length (149 um); it needs a length-preserving
repair or a check of the tracing/radius against the image.

## Length-preserving rod repair (2026-11-24): `rod_2026-11-24/`

`mesh_generation/tools/repair_trace_rod.py` (inextensible rod, bend 32 um,
clearance 55 um) keeps 97-98 % of each trace's length instead of contracting
it. `*_rod.swc/json/png` are single-trace repairs, `*_pad100_extended.swc` the
mesh inputs with frozen OuterWall tails, and `*_joint_with_*_coil.swc` the
joint third+fourth repair (55.1 um between the traces) for a shared Inter
mesh. The third-trace segment and full trace mesh r0-r2 with stable quality.
See `docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md`.
