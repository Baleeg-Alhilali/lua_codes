SOLVER CODE
===========

transport/ contains solute and membrane-transport models.
flow/ contains fluid and porous-flow models.
coupled/ contains multiphysics models.
archive/ contains superseded or experimental variants.

All active scripts calculate the Model repository root from their own location,
then resolve relative inputs and outputs from that root.

Example:
  /Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
    -ex code/flow/velocity_field.lua
