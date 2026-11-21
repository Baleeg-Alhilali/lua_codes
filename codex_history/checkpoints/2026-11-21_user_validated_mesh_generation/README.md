# User-validated mesh generation

This checkpoint preserves the complete mesh-generation source state selected by
the user on 2026-11-21. Generated meshes, solver results, logs, flow code,
transport code, and coupled code are intentionally excluded.

## Default single-nephron configuration

- SWC: `Nephron traces/3rd_trace_243_third.swc`
- lumen radius: 10 um
- membrane thickness: 10 um
- centerline scale: 24
- resample spacing: 5 um
- terminal spacing: 20 um
- terminal padding factor: 0
- anisotropy: 4
- bend anisotropy: 1
- bend chord-error factor: 0.01
- padding: 200 um
- BoxRefs: 4
- O-grid vertices: 10
- Inter topology: prism
- near-Inter prism thickness: 20 um
- Inter quality: 2
- refinement levels written: 2

The hybrid Inter topology uses prisms near the nephron, conforming
pyramid/tetrahedron transition cells, and tetrahedra in the far Inter domain.
Hexahedra remain in the structured Lumen/Membrane O-grid.

On branch `feature/inter-hexahedra`, `-interMesh hex` is also available. It
uses hexahedra in the structured near-Inter region and the same conforming
pyramid/tetrahedron transition. The default remains `prism`. See
`docs/INTER_HEXAHEDRA_VALIDATION_2026-11-21.md` for the R0--R2 validation and
the hex-specific refinement policy.

## Active source hashes

- `neurites_from_swc.cpp`: `e8f3b33f52a7f6556d15992977bf298f9e47b6c4153febf7976f4583bfe4df4b`
- `neurites_from_swc.h`: `0571c8966524fad903588fd71d1ef413d95f5ada4cc8549bed1e148ebb3ff60d`
- `neuro_collection_plugin.cpp`: `b91de5744b4306986e8dd2d0586997606d90871fd20449599051acff8f123fda`
- `swc_to_neurite_projector_nephron_final_pipeline.lua`: `e54a655280c41147fc3ae545990d366ac13e59f3765045f657cd7000b9e966cb`
- `quadratic_interpolate_swc.py`: `fac0831fcd8b50847bbd4f06412aa4e97faa91eaa744fc5e4406e5675e01426e`
- `extend_swc_terminals_to_outerwall.py`: `aaab2a5c3302f9ab664a9e0fcf575714bba8c6793559baf62c9a1948290b7477`
- `validate_ugx_mesh.py`: `a781042b0f533c6e75fd93e7b736eddb4f23ab8dd852fea0b8180669dcc03e62`
- `run_mesh_parameter_sweep.py`: `23937fc1d3ac44d4d6f221b3f44a66726d17935b4e4e053dc5171bdd9891a761`

The external active plugin source is not itself a Git worktree. Restore the
three `neuro_collection` files from this checkpoint into the active UG4 plugin
tree and rebuild `neuro_collection` to reproduce this source state.
