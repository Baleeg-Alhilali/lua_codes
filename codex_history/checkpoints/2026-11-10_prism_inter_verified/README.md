# Verified prismatic-Inter source checkpoint

This checkpoint preserves the exact modified `neuro_collection` sources used
for the 2026-11-10 straight-control r0--r2 validation. The active external
plugin directory is not a Git worktree, so this copy is the repository-backed
source record. The source and mesh hashes, command, validation counts, flow
diagnostics, and limitations are in
`docs/MESH_GENERATION_PRISM_INTER_VALIDATION_2026-11-10.md`.

The local checkpoint's `logs/` directory contains the compact generation and
flow records; those generated logs are intentionally omitted from the ordinary
source commit. Large UGX and VTU outputs remain in ignored `runs/` and
`Results_final/` directories.
