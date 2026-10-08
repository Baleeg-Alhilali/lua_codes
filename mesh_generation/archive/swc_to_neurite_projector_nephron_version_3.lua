-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/mesh_generation/archive/.*$"),
    "Expected this archived script below Model/mesh_generation/archive: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Generate a conforming Lumen/Membrane/Inter nephron mesh with Gmsh.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

local script_dir = MODEL_ROOT
local cfg = {}
cfg.swc = util.GetParam("-swc", script_dir .. "/Nephron traces/3rd_trace_243_section_96_154.swc")
cfg.out = util.GetParam("-out", script_dir .. "/runs/third_trace_section_nephron_version_3.ugx")
cfg.lumenRadius = util.GetParamNumber("-lumenRadius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
cfg.padding = util.GetParamNumber("-padding", 10.0)
cfg.nearSize = util.GetParamNumber("-nearSize", 0.5)
cfg.farSize = util.GetParamNumber("-farSize", 12.0)

local function quote(value)
	return "'" .. tostring(value):gsub("'", "'\\''") .. "'"
end
local function strip_extension(path)
	return (path:gsub("%.[^/%.]+$", ""))
end
local function dirname(path)
	return path:match("^(.*)/[^/]*$") or "."
end
local function positive(name, value)
	assert(value > 0, name .. " must be > 0")
end

positive("lumenRadius", cfg.lumenRadius)
positive("membraneThickness", cfg.membraneThickness)
positive("padding", cfg.padding)
positive("nearSize", cfg.nearSize)
positive("farSize", cfg.farSize)
assert(cfg.farSize >= cfg.nearSize, "farSize must be >= nearSize")
os.execute("mkdir -p " .. quote(dirname(cfg.out)))

local msh = strip_extension(cfg.out) .. ".msh"
local generator = script_dir .. "/generate_nephron_version_3.py"
local command = table.concat({
	"python3", quote(generator),
	"--swc", quote(cfg.swc),
	"--msh", quote(msh),
	"--ugx", quote(cfg.out),
	"--lumen-radius", tostring(cfg.lumenRadius),
	"--membrane-thickness", tostring(cfg.membraneThickness),
	"--padding", tostring(cfg.padding),
	"--near-size", tostring(cfg.nearSize),
	"--far-size", tostring(cfg.farSize)
}, " ")

print("Version 3 conforming nephron build")
print("  SWC              = " .. cfg.swc)
print("  output UGX       = " .. cfg.out)
print("  lumen radius     = " .. cfg.lumenRadius)
print("  membrane         = " .. cfg.membraneThickness)
print("  box padding      = " .. cfg.padding)
print("  near/far size    = " .. cfg.nearSize .. " / " .. cfg.farSize)
local ok = os.execute(command)
assert(ok == true or ok == 0, "Gmsh generator failed")

local mesh = Mesh()
assert(LoadMesh(mesh, cfg.out), "UG4 validation load failed: " .. cfg.out)
local sh = mesh:subset_handler()
assert(sh:num_subsets() == 8, "Expected eight named subsets")
-- The converter has already performed the ordered equivalent of selecting and
-- closing the volumes first, followed by the boundary faces. Repeating that
-- operation after loading would overwrite and erase the boundary markers.
print("Done. Subsets: Lumen, Membrane, Inter, Apical, Basolateral, OuterWall, Inlet, Outlet.")
