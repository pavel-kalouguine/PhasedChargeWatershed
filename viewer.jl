# Command-line launcher for the interactive density viewer.
#
#   julia --project viewer.jl path/to/data.json [path/to/results.jld2]
#
# The optional second argument is a file written by prewatershed.jl; when it is given,
# the boundaries between the watershed basins are drawn as white lines over the density.
#
# Opens one window. The "Add a view" button clones the current settings into a new
# window. Windows can be closed independently; closing the last one ends the program.

using GLMakie, JLD2
import PhasedChargeWatershed: load_data, build_viewer, global_density_limits

if !(1 <= length(ARGS) <= 2)
    error("usage: julia --project viewer.jl <data.json> [results.jld2]")
end
pd = load_data(ARGS[1])

result = nothing
if length(ARGS) == 2
    isfile(ARGS[2]) || error("results file not found: $(ARGS[2])")
    result = load_object(ARGS[2])

    if result.phased_data.peaks != pd.peaks
        @warn "the results were computed from different phased data than $(ARGS[1]); " *
              "the basin boundaries will not match the density"
    end
end

climits = result === nothing ? global_density_limits(pd) : global_density_limits(result)

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
