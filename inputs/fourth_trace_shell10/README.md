# Full fourth trace, 10 um hex Inter shell of two 5 um layers (2026-11-26)

**Use `fourth_shell10_a15.ugx`** (anisotropy 1.5, better caps; see the last
section). `fourth_shell10.ugx` (anisotropy 2) is kept for reference.
Load the coarse projector-backed UGX and refine in memory.
Subset names identical to `inputs/fourth_trace_full_pad100/`. Geometry pass
r0-r2 only; no flow or transport solve yet.

**New geometry.** The trace was re-repaired so that the shell fits everywhere:
`repair_trace_rod.py "Nephron traces/4th_trace_167.swc" --min-curvature-radius 55
--clearance 75 --inter-shell 10 --ramp-steps 24 --lbfgs-iters 1000 --w-length 1e5`
(fixed ends; the default homotopy failed at R 51 um). Coil 1128 um (99.4 % of
raw 1134 um; the older relaxed trace was 942 um), min bend radius 55.2 um,
min non-local distance 76.9 um. With tails to the OuterWall (padding 100 um):
1437 um. The coil differs in shape from `fourth_trace_full_pad100`.
- Inlet centre (mm) 0.338568, 0.400919, 0.163024; Outlet centre (mm)
  0.600995, -0.010000, 0.204605 (y = -0.01 OuterWall plane). Lumen radius 10 um.

Mesh: `rod_repair_workflow.sh` with
`INTER_OPTS="-prismInterThickness 10 -interMinLayers 2 -interMaxLayers 2 -anisotropy 2"`,
i.e. radius 10, membrane 10, ogrid 12, BoxRefs 3, hex Inter 10 um / 2 layers,
straight-path anisotropy 2 (was 4). Shell applied at full 10 um at every
station (0 curvature-limited). 242 axial stations, 86,630 coarse volumes,
r1 698,824, r2 5,625,296.

Generator validation PASS r0, r1, r2. Checker mode 2: 0 inverted/zero volumes.

| min corner scaled J | r0 | r1 | r2 |
|---|---|---|---|
| tet | 0.0252 | 0.0061 | 0.0054 |
| hex | 0.962 | 0.988 | 0.996 |
| prism | 0.498 | 0.470 | 0.471 |
| pyramid | 0.266 | 0.024 | 0.022 |

Hex edge aspect ratio (longest / shortest, `mesh_generation/tools/hex_aspect_ratio.py`):

| Subset | r0 max | r1 max | r2 max | hexes > 10 at r2 |
|---|---|---|---|---|
| Inter | 6.55 | 7.23 | 8.21 | 0 |
| Membrane | 5.21 | 6.21 | 8.21 | 0 |
| Lumen | 10.07 | 11.89 | 15.68 | 144 of 185,088 |

Limits: lumen-core hexes next to the 2.6 um inner ring exceed 10 and grow
with refinement; the cap pyramids drop r0 -> r1 (0.27 -> 0.024) and then hold.

## `fourth_shell10_a15.ugx`: anisotropy 1.5 (recommended, 2026-11-26)

Same trace and shell; only `-anisotropy 1.5` instead of 2. 174 axial stations,
67,925 coarse volumes, r1 553,040, r2 4,449,232 (smaller than anisotropy 2).
Generator validation PASS r0-r2; checker mode 2, 0 inverted/zero volumes:

| min corner scaled J | r0 | r1 | r2 |
|---|---|---|---|
| tet | 0.0271 | 0.0221 | 0.0221 |
| hex | 0.962 | 0.989 | 0.997 |
| prism | 0.499 | 0.470 | 0.471 |
| pyramid | 0.326 | 0.109 | 0.101 |

Hex edge aspect ratio max r0 / r1 / r2: Inter 4.23 / 5.14 / 7.14, Membrane
4.15 / 5.14 / 7.14, Lumen 8.01 / 9.85 / 13.65 (96 of 132,864 lumen-core hexes
above 10 at r2; none in Inter or Membrane).

Cap/TetGen test (same trace, coarse variants checked r0-r2): smaller TetGen
cells near the shell (`-interNearEdgeFactor 1 -interNearDistanceFactor 3`)
gave 5-6x more volumes and ~9x more coarse tets below J 0.05 without
improving the caps; `-interQuality 1.4` changed nothing. Anisotropy 1.5 alone
raised the pyramid minimum 0.022 -> 0.101 and the tet minimum
0.0054 -> 0.0221 at r2. Runs: `runs/fourth_shell10_caps_2026-11-26/`.
