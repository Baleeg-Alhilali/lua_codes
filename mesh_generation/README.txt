MESH GENERATION
===============

Objective: convert SWC nephron traces into projector-backed UGX meshes, refine
existing meshes, and inspect mesh quality.

Final pipeline (only these files are kept here):
  swc_to_neurite_projector_nephron_final_pipeline.lua  single nephron
  swc_to_neurite_projector_multiple_nephrons.lua       several nephrons, one Inter box
  refine_neurite_projector.lua                         conforming refinement of a coarse UGX
  check_refinement_quality.lua, Mesh_quality.lua       quality checks
  quadratic_interpolate_swc.py, extend_swc_terminals_to_outerwall.py,
  validate_ugx_mesh.py                                 helpers called by the pipeline
  tools/                                               trace repair (repair_trace_rod.py,
                                                       rod_repair_*workflow.sh), trace audit
                                                       and plots, hex aspect ratio

meshes/<YYYY-MM-DD>_<short_description>/ holds the generated meshes that are
kept, each with a README. Load the coarse projector-backed UGX in the solver
and refine in memory. meshes/legacy_transport_meshes/ holds the inputs of the
older transport models in code/transport.

See COMMANDS.txt for command examples and parameter explanations.
