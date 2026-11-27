SOLVER CODE
===========

flow/       current Stokes-Darcy flow models (fourth trace FE/FV1, joint traces)
transport/  transport models (nephron diffusion check and the older 2D/3D models)
coupled/    flow-transport coupling (nephron_flow_transport_standalone.lua)

Only working scripts are kept here. Debugging and test scripts are deleted once
their findings are in a dated note in docs/; git history keeps old versions.
All scripts calculate the Model root from their own location. Meshes are read
from mesh_generation/meshes/; kept results are written to
Results/<YYYY-MM-DD>_<short_description>/.

Example:
  mpirun -np 8 /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
    -ex code/flow/nephron_flow_fourth_pad100.lua -numRefs 2
