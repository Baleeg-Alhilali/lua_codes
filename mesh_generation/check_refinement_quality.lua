ug_load_script("ug_util.lua")
InitUG(3, AlgebraType("CPU", 1))
local f = util.GetParam("-grid", "")
local n = util.GetParamNumber("-refs", 1)
local mode = util.GetParamNumber("-mode", 0)
check_local_nephron_refinement(f, n, mode,
  util.GetParamNumber("-cx", 0), util.GetParamNumber("-cy", 0), util.GetParamNumber("-cz", 0),
  util.GetParamNumber("-radius", 0), util.GetParamNumber("-buffer", 1))
print("CHECK DONE")
