# Model Lua scripts

This directory contains the active Lua scripts that were previously stored at
the `Model` directory root.

Each script resolves its location and changes the UG4 working directory back to
the parent `Model` directory before loading inputs or writing outputs. Existing
relative paths therefore continue to refer to folders such as:

- `grids/`
- `Nephron traces/`
- `ProMeshFiles/`
- `runs/` and `run/`
- `Results/`, `results/`, and `output/`

Scripts can be launched from the `Model` directory:

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell -ex code/velocity_field.lua
```

They can also be launched by absolute path from another directory:

```bash
/Users/alhilaba/UG4_promesh_ogrid/bin/ugshell \
  -ex /Users/alhilaba/Documents/Promesh_swc/Model/code/velocity_field.lua
```

Relative command-line paths such as `-grid runs/example.ugx` are interpreted
from the `Model` directory. Absolute command-line paths continue to work.

The historical Lua copies under `ProMeshFiles/` remain in place to avoid
overwriting active scripts with the same filenames.
