-- Resolve every relative input and output path from the Model directory.
local __model_script = debug.getinfo(1, "S").source
if __model_script:sub(1, 1) == "@" then __model_script = __model_script:sub(2) end
if __model_script:sub(1, 1) ~= "/" then
    __model_script = CurrentWorkingDirectory() .. "/" .. __model_script
end
local MODEL_ROOT = assert(__model_script:match("^(.*)/code/[^/]+$"),
    "Expected this script inside Model/code: " .. __model_script)
ChangeDirectory(MODEL_ROOT)
-- test for the github
-- test
-- Na+/K+-ATPase variant of the 2D membrane-transport model.
-- The shared implementation switches to Na/K fields and the compiled NaKPump
-- transporter when this flag is set before loading it.
naKPumpMode = true
combinedTransportMode = false
ug_load_script("code/2D_trans_combined_print.lua")
