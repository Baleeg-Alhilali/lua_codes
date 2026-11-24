# UG4 build overlay — state of 2026-11-24

This folder reproduces the exact UG4 build used for the true-scale fourth
trace work (mesh generator, NeuriteProjector, TetGen wrapper, pressure-trace
interface diagnostics). `UG4_promesh_ogrid` itself is not a git repository;
`ugcore` and the plugins are clones of the public UG4 repositories with local
changes, and `plugins/neuro_collection` is not under git at all. Everything
needed is therefore stored here.

## Contents

| Path | What |
|---|---|
| `patches/ugcore.patch` | all local changes of `ugcore` against `46b9b57a8` (15 files) |
| `patches/plugin_NavierStokes.patch` | local changes of `plugins/NavierStokes` against `8df0ce2` (5 files) |
| `neuro_collection/` | complete source of `plugins/neuro_collection` (no build products) |
| `externals_EigenForUG4_untracked/` | 9 files present in `externals/EigenForUG4` but untracked upstream |
| `ug4_root/CMakeLists.txt`, `ug4_root/dot_ughub/` | ug4 root files from the ughub setup |
| `SHA256SUMS_installed.txt` | hashes of every overlay-provided file as installed in the UG4 tree |
| `apply_overlay.sh` | applies everything to a UG4 tree, verifies hashes, optionally builds |

## Upstream commits

| Repository | Commit | Local changes |
|---|---|---|
| ugcore (github.com/UG4/ugcore) | 46b9b57a8c8398939bf1b6a4661396c863f6a723 | patch |
| plugin_NavierStokes | 8df0ce228312397f10638219ba39383dbf59b4ac | patch |
| plugin_ConvectionDiffusion | c55e90170fc1482ddf54f2989d3418b7dbb545a6 | none |
| plugin_Limex | cb076d34a8f7fb68aa139486216468c2adf31533 | none |
| plugin_LuaShell | f4bfa0e6e3a61459330817399161ef2b94ac2238 | none (disabled in build) |
| plugin_ProMesh | 9f3ee2ea2a4b27cab9be59e814bc57e354b1a8e3 | none |
| plugin_SuperLU6 | 50b43fdac91aba42057afedafdfb0d32e6c099ce | none |
| external_BoostForUG4 | dbdd3c8cc7079c056f022d5f22545f18d5edf0f5 | none |
| external_EigenForUG4 | 0ebff70a3d12b664beefab1b5b9b73515330a288 | 9 untracked files (copied) |
| neuro_collection | not git | full copy |

Key post-change hashes (also in `SHA256SUMS_installed.txt`):

- `ugcore/.../projectors/neurite_projector.h` `914862019aba9df74b101018df5de015f0f6caa7c426931dd2c738c6b72ee280`
- `ugcore/.../projectors/neurite_projector.cpp` `5370a773ff6264b4aa03f312848bbf54ec023ad65dac829cf34a699fc60fe941`
- `ugcore/.../grid_generation/tetrahedralization.cpp` `51b0b40b4436155fe3b1c063bba7ba612e9dd0c963d079afecce3585721e2d62`
- `neuro_collection/grid_generation/neurites_from_swc.cpp` `55402eebee844c48417a62f237b8db3c3a2d4c6f7ff5aac8e15f172aff5b3fab` (2026-11-24 linear outward caps + checker lineage/mode 4; previous `c61a0c30…`)
- `neuro_collection/pressure_trace_constraint.cpp` `5e0d776fb4e939036ebebb21f2094f3e82ea99af9f686354a040526f04b8e63c`

## Build configuration used (this workstation)

macOS 27.2.1, Apple clang 21.0.0, CMake 4.3.2, Open MPI 5.0.9.
CMake: Release, DIM=ALL, CPU=ALL, PRECISION=double, PARALLEL=ON, LAPACK=ON,
BLAS=ON (Accelerate), INTERNAL_BOOST=ON, LINK_TETGEN=ON, tetgen=ON,
COMPILE_INFO=ON, STATIC_BUILD=OFF, EMBEDDED_PLUGINS=OFF; plugins ON:
ConvectionDiffusion, Limex, NavierStokes, ProMesh, SuperLU6
(USE_INTERNAL_SUPERLU=ON), neuro_collection; LuaShell OFF.

## Runtime notes

- The nephron importer sets `UG_TETGEN_DOUBLE_PRECISION=1` unless the
  variable is already set (core default stays single precision).
- `UG_TETGEN_MIN_DIHEDRAL` (degrees) is an opt-in TetGen quality bound; not
  used for the reported meshes.
- NeuriteProjector serializes a radial-fade pair from class version 1; UGX
  files written with this build need this build (older files still load).
- `UG_NEPHRON_PROJECTED_CAPS=1` restores the fully projected outward caps
  (default since 2026-11-24: caps refine linearly, base quads stay on the tube
  map; see `docs/ROD_REPAIR_AND_LINEAR_CAPS_2026-11-24.md`).
