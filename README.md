# Nephron model code

This repository contains the Lua, Python, documentation, and compact input
files used for the nephron simulations in UG4. The local repository root is:

`/Users/alhilaba/Documents/Promesh_swc/Model`

The repository is intentionally separated from generated meshes and numerical
results. Kept results are in `Results/<YYYY-MM-DD>_<desc>/` (ignored by Git); kept
meshes are in `mesh_generation/meshes/<YYYY-MM-DD>_<desc>/`.

Primary source areas:

- `code/transport/` — solute transport and electro-diffusion models.
- `code/flow/` — Navier–Stokes, Darcy, and velocity-field models.
- `code/coupled/` — coupled flow/transport and multicomponent models.
- `mesh_generation/` — the final mesh pipeline, its tools and `meshes/`.
- `docs/` — experiment notes, handoffs, and project documentation.
- `Nephron traces/` and `grids/` — compact source inputs.

Start with `FOLDER_STRUCTURE.txt`, `BRANCHES.txt`, and
`mesh_generation/COMMANDS.txt`.
