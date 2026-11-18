# Nephron model code

This repository contains the Lua, Python, documentation, and compact input
files used for the nephron simulations in UG4. The local repository root is:

`/Users/alhilaba/Documents/Promesh_swc/Model`

The repository is intentionally separated from generated meshes and numerical
results. Large or reproducible outputs remain on the workstation under
`runs/`, `Results/`, `results/`, and `output/`; those directories are ignored
by Git.

Primary source areas:

- `code/transport/` — solute transport and electro-diffusion models.
- `code/flow/` — Navier–Stokes, Darcy, and velocity-field models.
- `code/coupled/` — coupled flow/transport and multicomponent models.
- `code/archive/` — older or experimental solver variants retained for study.
- `mesh_generation/` — active mesh-generation pipelines and helper tools.
- `docs/` — experiment notes, handoffs, and project documentation.
- `Nephron traces/` and `grids/` — compact source inputs.

Start with `FOLDER_STRUCTURE.txt`, `BRANCHES.txt`, and
`mesh_generation/COMMANDS.txt`.
