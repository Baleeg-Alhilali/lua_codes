# Full fourth trace, true scale, 100 um padding (2026-11-24)

- `fourth_trace_full_pad100.ugx`: coarse projector-backed mesh used for all
  2026-11-24 flow results (load this; refine in memory). Regenerating it is
  possible but TetGen output varies slightly between runs.
- `fourth_trace_full_repaired_centreline.swc`: repaired coil (um), from
  `Nephron traces/4th_trace_167.swc`, 1 um voxels:
  `python3 mesh_generation/tools/repair_trace_clearance.py "Nephron traces/4th_trace_167.swc" --clearance 55 -o <out>`
- `fourth_trace_full_pad100_extended_repaired.swc`: with OuterWall terminal
  extensions (padding 100 um) and tails frozen 40 um, repaired again.
- Regenerate the mesh from Model/mesh_generation:
  `<UG4>/bin/ugshell -ex swc_to_neurite_projector_nephron_final_pipeline.lua -- -swc ../inputs/fourth_trace_full_pad100/fourth_trace_full_pad100_extended_repaired.swc -skipSmoothing 1 -centerlineScale 1 -extendTerminals 0 -radius 10 -membraneThickness 10 -ogrid 12 -BoxRefs 3 -interMesh hex -prismInterThickness 3 -numRefs 2 -out <out>.ugx`
- Inlet centre (mm) 0.360781065, 0.393057145, 0.162855950, direction +x;
  lumen radius 10 um; lumen path 1.206 mm.
