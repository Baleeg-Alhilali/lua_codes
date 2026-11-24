# Full third trace, true scale, 100 um padding (2026-11-24)

- `third_trace_full_pad100.ugx`: coarse projector-backed mesh (load this and
  refine in memory). Generated with linear outward caps (neurites_from_swc.cpp
  `55402eeb...`); r0-r2 generator validation PASS; checker mode 2: tet 0.0266 /
  0.0266, hex 0.966 / 0.978, pyramid 0.146 / 0.136 at r1 / r2. 23,982 vertices,
  73,379 volumes. Geometry pass only; no flow solve yet.
- `third_trace_full_repaired_centreline.swc`: length-preserving rod repair of
  `Nephron traces/3rd_trace_243.swc` (1318.8 -> 1293.5 um, bend 32, clearance 55).
- `third_trace_full_pad100_extended_repaired.swc`: OuterWall terminal
  extensions (padding 100 um), frozen 40 um tails, rod-repaired again
  (1585 um lumen path).
- Regenerate (from Model):
  `mesh_generation/tools/rod_repair_workflow.sh "Nephron traces/3rd_trace_243.swc" <out_dir> <prefix>`
  (TetGen output varies slightly between runs).
- See `docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md`.
