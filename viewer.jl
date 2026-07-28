# Command-line launcher for the interactive density viewer.
#
#   julia --project viewer.jl path/to/data.json
#   julia --project viewer.jl path/to/results.jld2
#
# What is drawn follows from the file that is given. A .json holds phased data only, so
# the density alone is shown; a .jld2 written by prewatershed.jl holds a watershed result
# as well, and the boundaries between the basins are drawn as white lines over the
# density, with a row of controls to postprocess them. A result carries the phased data it
# was computed from, so the density and the basins on screen always belong together.
#
# Opens one window. The "Add a view" button clones the current settings into a new
# window. Windows can be closed independently; closing the last one ends the program.

using GLMakie
import PhasedChargeWatershed: load_data, load_result, build_viewer, global_density_limits

length(ARGS) == 1 || error("usage: julia --project viewer.jl <data.json | results.jld2>")
input_path = ARGS[1]
isfile(input_path) || error("input file not found: $input_path")

extension = lowercase(splitext(input_path)[2])
if extension == ".json"
    pd = load_data(input_path)
    result = nothing
    climits = global_density_limits(pd)
elseif extension == ".jld2"
    watershed = load_result(input_path)
    pd = watershed.phased_data
    climits = global_density_limits(watershed)
    result = Observable(watershed)
else
    error("do not know what to do with \"$extension\": " *
          "expected .json (phased data) or .jld2 (watershed results)")
end

# The set of open windows; the program runs until the last one is closed
screens = Set{GLMakie.Screen}()

function add_view(init = nothing)
    fig = build_viewer(pd; on_add_view = add_view, init = init, colorrange = climits,
                       result = result)
    screen = GLMakie.Screen()
    display(screen, fig)                 # open a new window
    push!(screens, screen)
    on(events(fig.scene).window_open) do isopen
        isopen || delete!(screens, screen)   # drop it once the window closes
    end
    return
end

add_view()
while !isempty(screens)   # running until every window is closed
    sleep(0.1)
end
