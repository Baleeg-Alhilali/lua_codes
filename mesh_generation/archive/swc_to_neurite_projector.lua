-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/mesh_generation/archive/.*$"),
    "Expected this archived script below Model/mesh_generation/archive: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Build a projector-backed neurite surface directly from an SWC tree.
--
-- This deliberately does not embed the morphology in a box.  It uses UG4's
-- existing neuro_collection importer, which writes a .ugx containing a
-- ProjectionHandler with a serialized NeuriteProjector.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

local cfg = {}
cfg.swc = util.GetParam("-swc", "sample_tree.swc")
cfg.out = util.GetParam("-out", "runs/swc_projector.ugx")
cfg.radius = util.GetParamNumber("-radius", -1.0) -- < 0 keeps SWC radii
cfg.anisotropy = util.GetParamNumber("-anisotropy", 2.0)
cfg.numRefs = math.floor(util.GetParamNumber("-numRefs", 1))

local function dirname(path)
	local dir = path:match("^(.*)/[^/]*$")
	if dir == nil or dir == "" then
		return "."
	end
	return dir
end

local function strip_extension(path)
	return (path:gsub("%.[^/%.]+$", ""))
end

local function shell_quote(s)
	return "'" .. tostring(s):gsub("'", "'\\''") .. "'"
end

local function ensure_dir(path)
	if path ~= "." then
		local ok = os.execute("mkdir -p " .. shell_quote(path))
		assert(ok == true or ok == 0, "Could not create directory: " .. path)
	end
end

local function assert_positive(name, value)
	assert(value > 0.0, name .. " must be > 0")
end

local function radius_tag(radius)
	local tag = string.format("%.8g", radius)
	tag = tag:gsub("%.", "p"):gsub("%-", "m"):gsub("%+", "")
	return tag
end

local function write_radius_override_swc(src, dst, radius)
	local in_file = assert(io.open(src, "r"), "Cannot open SWC input: " .. src)
	local out_file = assert(io.open(dst, "w"), "Cannot write SWC radius override: " .. dst)

	for line in in_file:lines() do
		local hash_pos = line:find("#", 1, true)
		local payload = line
		local comment = ""
		if hash_pos then
			payload = line:sub(1, hash_pos - 1)
			comment = line:sub(hash_pos)
		end

		local tokens = {}
		for tok in payload:gmatch("%S+") do
			tokens[#tokens + 1] = tok
		end

		if #tokens == 0 then
			out_file:write(line, "\n")
		else
			assert(#tokens == 7, "SWC line must have 7 columns: " .. line)
			local out_line = table.concat({
				tokens[1],
				tokens[2],
				tokens[3],
				tokens[4],
				tokens[5],
				string.format("%.17g", radius),
				tokens[7]
			}, " ")
			if comment ~= "" then
				out_line = out_line .. " " .. comment
			end
			out_file:write(out_line, "\n")
		end
	end

	in_file:close()
	out_file:close()
end

local function file_contains(path, needle)
	local file = assert(io.open(path, "r"), "Cannot open file: " .. path)
	local text = file:read("*a")
	file:close()
	return text:find(needle, 1, true) ~= nil
end

local function assert_neurite_projector_ugx(path)
	assert(file_contains(path, "<projector type=\"neurite\""),
		"UGX is missing a serialized NeuriteProjector: " .. path)
	assert(file_contains(path, "<vertex_attachment name=\"npSurfParams\""),
		"UGX is missing npSurfParams. Refinement cannot project neurite vertices: " .. path)
end

local function refine_with_promesh(src, dst, refs)
	local mesh = Mesh()
	assert(LoadMesh(mesh, src), "LoadMesh failed: " .. src)
	for i = 1, refs do
		SelectAll(mesh)
		print("  projected refine " .. i .. " / " .. refs)
		Refine(mesh)
	end
	assert(SaveMesh(mesh, dst), "SaveMesh failed: " .. dst)
end

assert(cfg.swc ~= "", "Pass an SWC file with -swc <path>")
assert(cfg.out ~= "", "Pass an output path with -out <path>")
assert_positive("anisotropy", cfg.anisotropy)
assert(cfg.numRefs >= 0, "numRefs must be >= 0")
if cfg.radius == 0.0 then
	error("radius must be > 0, or < 0 to keep the SWC radii")
end

ensure_dir(dirname(cfg.out))

local import_swc = cfg.swc
if cfg.radius > 0.0 then
	assert_positive("radius", cfg.radius)
	import_swc = strip_extension(cfg.out) .. "_radius_" .. radius_tag(cfg.radius) .. ".swc"
	write_radius_override_swc(cfg.swc, import_swc, cfg.radius)
end

local out_ugx = strip_extension(cfg.out) .. ".ugx"
local coarse_ugx = out_ugx
if cfg.numRefs > 0 then
	coarse_ugx = strip_extension(cfg.out) .. "_coarse.ugx"
end

print("SWC projector build")
print("  input SWC  = " .. cfg.swc)
print("  import SWC = " .. import_swc)
print("  output UGX = " .. out_ugx)
if coarse_ugx ~= out_ugx then
	print("  coarse UGX = " .. coarse_ugx)
end
print("  radius     = " .. (cfg.radius > 0.0 and tostring(cfg.radius) or "from SWC"))
print("  anisotropy = " .. cfg.anisotropy)
print("  numRefs    = " .. cfg.numRefs)

-- Keep importer refinement disabled here.  Its refined hierarchy files are
-- useful for solver hierarchies, but they do not carry defPH in the final UGX.
import_neurites_from_swc(import_swc, coarse_ugx, cfg.anisotropy, 0)
assert_neurite_projector_ugx(coarse_ugx)

if cfg.numRefs > 0 then
	refine_with_promesh(coarse_ugx, out_ugx, cfg.numRefs)
	assert_neurite_projector_ugx(out_ugx)
end

print("Done. The UGX contains defPH / NeuriteProjector / npSurfParams data for refinement.")
