-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/mesh_generation/archive/.*$"),
    "Expected this archived script below Model/mesh_generation/archive: " .. __model_script)
ChangeDirectory(MODEL_ROOT)

-- Build a two-domain nephron mesh directly from the specified SWC tree.
--
-- Lumen is the inner volume, Membrane is the surrounding hollow volume,
-- Apical is their shared interface, and Basolateral is the outer surface.

ug_load_script("ug_util.lua")
AssertPluginsLoaded({"neuro_collection", "ProMesh"})
InitUG(3, AlgebraType("CPU", 1))

local script_dir = MODEL_ROOT

local cfg = {}
cfg.swc = util.GetParam("-swc", script_dir .. "/Nephron traces/3rd_trace_243_section_96_154.swc")
cfg.out = util.GetParam("-out", script_dir .. "/runs/third_trace_section_nephron_version_2.ugx")
cfg.lumenRadius = util.GetParamNumber("-lumenRadius", 1.0)
cfg.membraneThickness = util.GetParamNumber("-membraneThickness", 0.2)
cfg.anisotropy = util.GetParamNumber("-anisotropy", 2.0)
cfg.numRefs = math.floor(util.GetParamNumber("-numRefs", 1))

-- Surface locations where a tube twist was observed. These diagnostics do not
-- modify the SWC path or the generated geometry.
cfg.twistProbes = {
	{318.979, 186.703, 57.9855},
	{318.143, 192.177, 59.916}
}
cfg.twistSmoothHalfWindow = 30
cfg.twistSmoothPasses = 60
cfg.twistSmoothStrength = 0.5

-- Replace only the known tight zigzag in the generated SWC. Point IDs outside
-- this interval retain their original coordinates exactly.
cfg.problemSegments = {
}

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

local vec_sub
local vec_length

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
				(tokens[7] == "-1") and "1" or (tokens[2] == "0" and "3" or tokens[2]),
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

local function replace_problem_segments(path, segments)
	local points = {}
	local by_id = {}
	local input = assert(io.open(path, "r"), "Cannot open generated SWC: " .. path)
	local lines = {}
	for line in input:lines() do
		local entry = {original = line}
		local hash_pos = line:find("#", 1, true)
		local payload = hash_pos and line:sub(1, hash_pos - 1) or line
		entry.comment = hash_pos and line:sub(hash_pos) or ""
		local tokens = {}
		for token in payload:gmatch("%S+") do tokens[#tokens + 1] = token end
		if #tokens > 0 then
			assert(#tokens == 7, "SWC line must have 7 columns: " .. line)
			entry.tokens = tokens
			entry.id = tonumber(tokens[1])
			entry.pos = {tonumber(tokens[3]), tonumber(tokens[4]), tonumber(tokens[5])}
			entry.pointIndex = #points + 1
			points[#points + 1] = entry
			by_id[entry.id] = entry
		end
		lines[#lines + 1] = entry
	end
	input:close()

	local function unit_direction(a, b)
		local direction = vec_sub(b, a)
		local length = vec_length(direction)
		assert(length > 1e-12, "Cannot construct a tangent from duplicate SWC points")
		return {direction[1] / length, direction[2] / length, direction[3] / length}
	end

	for _, segment in ipairs(segments) do
		local first = assert(by_id[segment.startId],
			"Problem-segment start ID not found: " .. segment.startId)
		local last = assert(by_id[segment.endId],
			"Problem-segment end ID not found: " .. segment.endId)
		local first_ind, last_ind = first.pointIndex, last.pointIndex
		assert(first_ind > 1 and last_ind < #points and first_ind < last_ind,
			"Problem segment must have unchanged neighboring points")
		assert(segment.handleScale > 0.0, "Bezier handleScale must be > 0")

		local start_tangent = unit_direction(points[first_ind - 1].pos, first.pos)
		local end_tangent = unit_direction(last.pos, points[last_ind + 1].pos)
		local chord = vec_length(vec_sub(last.pos, first.pos))
		local handle_length = segment.handleScale * chord
		local control1, control2 = {}, {}
		for d = 1, 3 do
			control1[d] = first.pos[d] + handle_length * start_tangent[d]
			control2[d] = last.pos[d] - handle_length * end_tangent[d]
		end

		local max_displacement = 0.0
		for i = first_ind + 1, last_ind - 1 do
			local t = (i - first_ind) / (last_ind - first_ind)
			local s = 1.0 - t
			local old = points[i].pos
			local corrected = {}
			for d = 1, 3 do
				corrected[d] = s * s * s * first.pos[d]
					+ 3.0 * s * s * t * control1[d]
					+ 3.0 * s * t * t * control2[d]
					+ t * t * t * last.pos[d]
			end
			max_displacement = math.max(max_displacement,
				vec_length(vec_sub(corrected, old)))
			points[i].pos = corrected
		end
		print(string.format(
			"  corrected SWC IDs %d..%d with cubic segment; max displacement %.6g",
			segment.startId, segment.endId, max_displacement))
	end

	local output = assert(io.open(path, "w"), "Cannot write corrected SWC: " .. path)
	for _, entry in ipairs(lines) do
		if entry.tokens then
			entry.tokens[3] = string.format("%.17g", entry.pos[1])
			entry.tokens[4] = string.format("%.17g", entry.pos[2])
			entry.tokens[5] = string.format("%.17g", entry.pos[3])
			local line = table.concat(entry.tokens, " ")
			if entry.comment ~= "" then line = line .. " " .. entry.comment end
			output:write(line, "\n")
		else
			output:write(entry.original, "\n")
		end
	end
	output:close()
end

vec_sub = function(a, b)
	return {a[1] - b[1], a[2] - b[2], a[3] - b[3]}
end

local function vec_dot(a, b)
	return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
end

local function vec_cross(a, b)
	return {
		a[2] * b[3] - a[3] * b[2],
		a[3] * b[1] - a[1] * b[3],
		a[1] * b[2] - a[2] * b[1]
	}
end


vec_length = function(v)
	return math.sqrt(vec_dot(v, v))
end

local function angle_deg(a, b)
	local denom = vec_length(a) * vec_length(b)
	if denom < 1e-12 then return nil end
	local c = vec_dot(a, b) / denom
	c = math.max(-1.0, math.min(1.0, c))
	return math.deg(math.acos(c))
end

local function read_swc_points(path)
	local points = {}
	local file = assert(io.open(path, "r"), "Cannot open SWC for twist diagnostics: " .. path)
	for line in file:lines() do
		if not line:match("^%s*#") then
			local tokens = {}
			for tok in line:gmatch("%S+") do tokens[#tokens + 1] = tok end
			if #tokens == 7 then
				points[#points + 1] = {
					id = tonumber(tokens[1]),
					pos = {tonumber(tokens[3]), tonumber(tokens[4]), tonumber(tokens[5])}
				}
			end
		end
	end
	file:close()
	return points
end

local function smooth_twist_probe_regions(path, probes, half_window, passes, strength)
	local points = read_swc_points(path)
	assert(#points >= 4, "Cannot smooth an SWC with fewer than four points")
	assert(half_window >= 1, "twist smoothing half-window must be >= 1")
	assert(passes >= 1, "twist smoothing passes must be >= 1")
	assert(strength > 0.0 and strength <= 1.0,
		"twist smoothing strength must be in (0, 1]")

	local selected = {}
	for _, probe in ipairs(probes) do
		local nearest_ind, nearest_dist = 1, math.huge
		for i, point in ipairs(points) do
			local dist = vec_length(vec_sub(point.pos, probe))
			if dist < nearest_dist then
				nearest_ind, nearest_dist = i, dist
			end
		end

		local first = math.max(2, nearest_ind - half_window)
		local last = math.min(#points - 1, nearest_ind + half_window)
		for i = first, last do selected[i] = true end
		print(string.format(
			"  smoothing around probe -> nearest SWC id %d, indices %d..%d",
			points[nearest_ind].id, first, last))
	end

	local original = {}
	for i, point in ipairs(points) do
		original[i] = {point.pos[1], point.pos[2], point.pos[3]}
	end

	for _ = 1, passes do
		local next_pos = {}
		for i, point in ipairs(points) do
			next_pos[i] = {point.pos[1], point.pos[2], point.pos[3]}
		end
		for i in pairs(selected) do
			for d = 1, 3 do
				local neighbor_avg = 0.5 * (points[i - 1].pos[d] + points[i + 1].pos[d])
				next_pos[i][d] = (1.0 - strength) * points[i].pos[d]
					+ strength * neighbor_avg
			end
		end
		for i, point in ipairs(points) do point.pos = next_pos[i] end
	end

	local max_displacement = 0.0
	for i in pairs(selected) do
		max_displacement = math.max(max_displacement,
			vec_length(vec_sub(points[i].pos, original[i])))
	end

	local input = assert(io.open(path, "r"), "Cannot reopen SWC after smoothing: " .. path)
	local lines = {}
	for line in input:lines() do lines[#lines + 1] = line end
	input:close()

	local point_ind = 0
	for line_ind, line in ipairs(lines) do
		local hash_pos = line:find("#", 1, true)
		local payload = hash_pos and line:sub(1, hash_pos - 1) or line
		local comment = hash_pos and line:sub(hash_pos) or ""
		local tokens = {}
		for tok in payload:gmatch("%S+") do tokens[#tokens + 1] = tok end
		if #tokens == 7 then
			point_ind = point_ind + 1
			local pos = points[point_ind].pos
			tokens[3] = string.format("%.17g", pos[1])
			tokens[4] = string.format("%.17g", pos[2])
			tokens[5] = string.format("%.17g", pos[3])
			lines[line_ind] = table.concat(tokens, " ")
			if comment ~= "" then lines[line_ind] = lines[line_ind] .. " " .. comment end
		end
	end

	local output = assert(io.open(path, "w"), "Cannot write smoothed SWC: " .. path)
	for _, line in ipairs(lines) do output:write(line, "\n") end
	output:close()
	print(string.format(
		"  localized twist smoothing: %d passes, max displacement %.6g",
		passes, max_displacement))
end

local function report_twist_probe_regions(path, probes)
	local points = read_swc_points(path)
	if #points < 4 then
		print("Twist diagnostics skipped: fewer than four SWC points.")
		return
	end

	local overall = vec_sub(points[#points].pos, points[1].pos)
	local ref_axis = 1
	if math.abs(overall[2]) < math.abs(overall[ref_axis]) then ref_axis = 2 end
	if math.abs(overall[3]) < math.abs(overall[ref_axis]) then ref_axis = 3 end
	local axis_names = {"X", "Y", "Z"}

	print("Twist-risk diagnostics")
	print("  projector reference axis = " .. axis_names[ref_axis])
	for probe_ind, probe in ipairs(probes) do
		local nearest_ind, nearest_dist = 1, math.huge
		for i, point in ipairs(points) do
			local dist = vec_length(vec_sub(point.pos, probe))
			if dist < nearest_dist then
				nearest_ind, nearest_dist = i, dist
			end
		end

		local first = math.max(3, nearest_ind - 10)
		local last = math.min(#points - 1, nearest_ind + 10)
		local max_turn, max_plane_rotation, min_segment = 0.0, 0.0, math.huge
		local max_ref_alignment = 0.0
		for i = first, last do
			local prev_tangent = vec_sub(points[i - 1].pos, points[i - 2].pos)
			local tangent = vec_sub(points[i].pos, points[i - 1].pos)
			local next_tangent = vec_sub(points[i + 1].pos, points[i].pos)
			min_segment = math.min(min_segment, vec_length(tangent), vec_length(next_tangent))
			local turn = angle_deg(tangent, next_tangent)
			if turn then max_turn = math.max(max_turn, turn) end
			local tangent_length = vec_length(tangent)
			if tangent_length > 1e-12 then
				max_ref_alignment = math.max(max_ref_alignment,
					math.abs(tangent[ref_axis]) / tangent_length)
			end

			local prev_normal = vec_cross(prev_tangent, tangent)
			local next_normal = vec_cross(tangent, next_tangent)
			if vec_length(prev_normal) > 1e-8 and vec_length(next_normal) > 1e-8 then
				local plane_rotation = angle_deg(prev_normal, next_normal)
				if plane_rotation then
					max_plane_rotation = math.max(max_plane_rotation, plane_rotation)
				end
			end
		end

		local point = points[nearest_ind]
		local min_ref_angle = math.deg(math.acos(math.min(1.0, max_ref_alignment)))
		local risk = max_ref_alignment > 0.95
		print(string.format(
			"  probe %d -> SWC id %d center (%.6g, %.6g, %.6g), surface distance %.6g",
			probe_ind, point.id, point.pos[1], point.pos[2], point.pos[3], nearest_dist))
		print(string.format(
			"    min reference angle %.2f deg, max turn %.2f deg, plane rotation %.2f deg, min segment %.6g -> %s",
			min_ref_angle, max_turn, max_plane_rotation, min_segment,
			risk and "TWIST RISK" or "reference-vector clearance OK"))
	end
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
assert_positive("lumenRadius", cfg.lumenRadius)
assert_positive("membraneThickness", cfg.membraneThickness)

local outer_radius = cfg.lumenRadius + cfg.membraneThickness
local lumen_scale = cfg.lumenRadius / outer_radius
ensure_dir(dirname(cfg.out))

local import_swc = cfg.swc
import_swc = strip_extension(cfg.out) .. "_outer_radius_" .. radius_tag(outer_radius) .. ".swc"
write_radius_override_swc(cfg.swc, import_swc, outer_radius)
replace_problem_segments(import_swc, cfg.problemSegments)

local out_ugx = strip_extension(cfg.out) .. ".ugx"
local coarse_ugx = out_ugx
if cfg.numRefs > 0 then
	coarse_ugx = strip_extension(cfg.out) .. "_coarse.ugx"
end

print("Two-domain nephron projector build")
print("  input SWC  = " .. cfg.swc)
print("  import SWC = " .. import_swc)
print("  output UGX = " .. out_ugx)
if coarse_ugx ~= out_ugx then
	print("  coarse UGX = " .. coarse_ugx)
end
print("  lumen radius       = " .. cfg.lumenRadius)
print("  membrane thickness = " .. cfg.membraneThickness)
print("  outer radius       = " .. outer_radius)
print("  lumen scale        = " .. lumen_scale)
print("  anisotropy = " .. cfg.anisotropy)
print("  numRefs    = " .. cfg.numRefs)

report_twist_probe_regions(import_swc, cfg.twistProbes)

-- Keep importer refinement disabled here. Its refined hierarchy files are
-- useful for solver hierarchies, but they do not carry defPH in the final UGX.
import_er_neurites_from_swc(import_swc, coarse_ugx,
	lumen_scale, cfg.anisotropy, 0)
assert_neurite_projector_ugx(coarse_ugx)

-- The ER importer already creates the required conforming core/shell topology.
-- Rename its four subsets to nephron terminology.
local subset_mesh = Mesh()
assert(LoadMesh(subset_mesh, coarse_ugx), "LoadMesh failed: " .. coarse_ugx)
SetSubsetName(subset_mesh, 0, "Membrane")
SetSubsetName(subset_mesh, 1, "Lumen")
SetSubsetName(subset_mesh, 2, "Basolateral")
SetSubsetName(subset_mesh, 3, "Apical")
assert(SaveMesh(subset_mesh, coarse_ugx), "SaveMesh failed: " .. coarse_ugx)
assert_neurite_projector_ugx(coarse_ugx)

if cfg.numRefs > 0 then
	refine_with_promesh(coarse_ugx, out_ugx, cfg.numRefs)
	assert_neurite_projector_ugx(out_ugx)
end

print("Done. Subsets: Lumen, Membrane, Apical, Basolateral.")
