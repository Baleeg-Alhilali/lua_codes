# New-build lumen-flow validation — 2026-10-05

## Files

- Whole generated mesh: `runs/3rd_trace_new_build_test.ugx`
- Regular lumen solver mesh: `runs/3rd_trace_new_build_test_lumen.ugx`
- Lumen extraction script: `extract_lumen_new_build.lua`
- Flow script: `flow_projector_new_build.lua`

The whole mesh is preserved unchanged. The lumen companion is produced in Lua
by removing the Membrane and Inter volumes (and their now-unused constrained
objects). It retains the serialized NeuriteProjector instances and
`npSurfParams`, so refinements are generated in memory from the coarse mesh.

No UG4 core or NavierStokes plugin source is modified. The flow setup is entirely
in Lua; the mesh geometry changes remain confined to neuro_collection.

## Proven solver configuration

The configuration follows the previously successful projector-flow logs:

- stationary Stokes FV1 in `Lumen`;
- Poiseuille inflow on `Inlet`, no-normal-stress `Outlet`, and no-slip
  `Apical,TripleJunction`;
- geometric multigrid using the actual level hierarchy;
- damped ILU (`0.7`) on refined levels;
- rediscretized level operators (`RAP=false`) and `StdTransfer`;
- gathered SuperLU only on level 0;
- all requested projector refinements are created before MPI distribution.

## Validation results

| Test | DoFs at top level | Linear iterations | Solve time | Result |
|---|---:|---:|---:|---|
| serial, ref1 | 321,276 | 6 | 5.926 s | PASS |
| serial, ref2 | 2,360,700 | 8 | 67.106 s | PASS |
| MPI 4 ranks, ref1 | 321,204 | 13 | 2.851 s | PASS |
| MPI 4 ranks, ref2 | 2,360,484 | 13 | 31.202 s | PASS |

Each run exited with status 0 and printed `FLOW_TEST_SUCCESS`.

## Output update — 2026-10-08

The validated solver configuration is unchanged. The Lua output setup was
updated so VTK output is enabled by default, files begin with `new_build_`, and
the results are written under:

`/Users/alhilaba/Documents/Promesh_swc/Model/Results`

The written fields are `velocity`, its separate `u`, `v`, and `w` components,
and `pressure`. Use `-writeVTK 0` only when an output-free test is intended.

## Reproduction

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex extract_lumen_new_build.lua -- \
  -grid runs/3rd_trace_new_build_test.ugx \
  -out runs/3rd_trace_new_build_test_lumen.ugx

/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex flow_projector_new_build.lua -- -numRefs 1

/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex flow_projector_new_build.lua -- -numRefs 2

mpirun -np 4 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex flow_projector_new_build.lua -- -numRefs 1

mpirun -np 4 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex flow_projector_new_build.lua -- -numRefs 2
```

## SHA-256

```text
877cb3f0841a4baf1246ec300e5a49052eca052578b9f4ed12a8bd5419dd84ef  flow_projector_new_build.lua
89b6806cdc5709d7fde96c50ea8c092ff2fc62a74920e62848adb423b051e0e2  extract_lumen_new_build.lua
b4fb1b1a7f3cf9d43eafd5a3a8397f66f4ab7109eb15d6d541a46d9f6875cb7f  runs/3rd_trace_new_build_test.ugx
bed27b29cff31559488aa22dfdcf06b8e8383dce7c543c70864433fd56f11487  runs/3rd_trace_new_build_test_lumen.ugx
```
