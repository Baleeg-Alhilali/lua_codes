# Codex project-history archive

Collected on 2026-11-07 for future Codex sessions. This directory contains
copies only; the original project files were not moved or changed.

## Session-history status

- `../sessions/` contains the copied VS Code Codex rollouts.
- No `archived_sessions` directory was present under this project or under
  `/Users/alhilaba/.codex` when this archive was assembled.
- `/Users/alhilaba/.codex/session_index.jsonl` has no archived-status field; it
  records only thread ID, thread name, and update time. Therefore none of the
  session files can be reliably classified as "archived" from that index.

## Archive layout

### `current_model/`

Copies from `/Users/alhilaba/Documents/Promesh_swc/Model`:

- `root_context/`: repository layout, branches, design decisions, and backup
  notes.
- `docs/`: velocity-field handoff, anisotropic-projector experiment history,
  and the flow-script structure contract.
- `mesh_generation/`: the projected-refinement defect handoff and its final
  validation results.
- `output/pdf/`: five technical reports covering flow, transport,
  Stokes-Darcy coupling, solver behavior, failed configurations, and handoff
  questions.
- `report_generators/`: the Python sources used to build those reports.
- `validation_logs_2026-10-13/`: raw level/refinement validation evidence.

### `legacy_promesh/`

Recovered from the earlier working root `/Users/alhilaba/Documents/ promesh`
(the space before `promesh` is part of the original path):

- the July nephron-plugin final reference;
- both O-grid working-build records;
- dated build, flow-validation, scale-adaptation, anisotropy, mesh-pipeline,
  and refinement-connectivity histories;
- the nephron+blood starter reference-build note; and
- the 20-page nephron geometry/refinement pipeline report.

The exact duplicate
`/Users/alhilaba/Documents/ promesh/tmp/MESH_REFINEMENT_VALIDATION_RESULTS_2026-10-13.md`
was not copied twice; its identical copy is in `current_model/mesh_generation/`.

### `legacy_onedrive/`

Contains the older OneDrive version of `VELOCITY_FIELD_HANDOFF.md`. The newer
`current_model/docs/VELOCITY_FIELD_HANDOFF.md` updates the machine paths and
UG4 command, so both versions are retained with their provenance separated.

## Recommended starting documents

For mesh-generation debugging, read in this order:

1. `legacy_promesh/work/MESH_GEOMETRY_PIPELINE_HISTORY_2026-10-08.md`
2. `legacy_promesh/work/outerwall_portal/REFINEMENT_CONNECTIVITY_FIX_2026-10-04.md`
3. `current_model/mesh_generation/MESH_REFINEMENT_VALIDATION_HANDOFF.md`
4. `current_model/mesh_generation/MESH_REFINEMENT_VALIDATION_RESULTS_2026-10-13.md`

For flow and transport work, start with:

1. `current_model/docs/VELOCITY_FIELD_HANDOFF.md`
2. `current_model/output/pdf/Urea_flow_projector_transport_technical_summary.pdf`
3. `current_model/output/pdf/working_flow_coupling_simple_long_tube_explained.pdf`
4. `current_model/output/pdf/stokes_darcy_coupling_solver_explanation.pdf`
5. `current_model/output/pdf/stokes_darcy_flow_coupling_handoff.pdf`

