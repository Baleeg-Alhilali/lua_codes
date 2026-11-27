# Third + fourth full traces in one shared Inter domain (2026-11-24)

- `third_fourth_joint_pad100.ugx`: coarse projector-backed two-nephron mesh
  (load this and refine in memory). Subsets `Lumen_i, Membrane_i, Apical_i,
  Basolateral_i, Inlet_i, Outlet_i, TripleJunction_i, MembraneWall_i` for
  i = 1 (third trace) and 2 (fourth trace); shared `Inter`, `OuterWall`.
  42,484 vertices, 122,420 volumes; hex Inter shell 3 um with linear outward
  caps; ogrid 12, BoxRefs 3, padding 100 um.
- Checks: coarse validation PASS; in-memory checker (mode 2, global) r1 / r2:
  tet 0.00222 / 0.00222, hex 0.971 / 0.979, prism 0.469 / 0.471, pyramid
  0.191 / 0.139, zero inverted cells (r2: 7.93 M leaf volumes, 11.7 GB peak).
  Generator r1 flat-file validation PASS on a sibling build of the same input
  (TetGen output differs slightly between runs). Geometry pass only: no flow or
  transport solve yet.
- Traces (um): joint length-preserving rod repair of `3rd_trace_243.swc` and
  `4th_trace_167.swc` (bend 35 um, clearance 55 um within and between traces,
  free ends): 1318.8 -> 1289.2 um and 1134.1 -> 1114.8 um.
  `third_fourth_joint_pad100_extended_repaired.swc` is the merged two-tree
  OuterWall-extended input (repaired again with 40 um frozen tails): lumen
  paths 1581 um (third) and 1507 um (fourth); 56.2 um between the traces.
- Terminals (mm): nephron 1 inlet on xmin at (0.17342, 0.43299, 0.23718),
  outlet on ymin at (0.60752, -0.01697, 0.18823); nephron 2 inlet on ymax at
  (0.52859, 0.63751, 0.18823), outlet on ymin at (0.60752, -0.01697, 0.28613).
- Regenerate (from Model):
  `mesh_generation/tools/rod_repair_joint_workflow.sh <out_dir> <prefix> "Nephron traces/3rd_trace_243.swc" "Nephron traces/4th_trace_167.swc"`
- See `docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md`, section 6.
