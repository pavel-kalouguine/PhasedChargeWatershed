# Command-line launcher for the interactive density viewer.
#
#   julia --project viewer.jl path/to/data.json
#
# Opens one window. The "Add a view" button clones the current settings into a new
# window. Windows can be closed independently; closing the last one ends the program.

using GLMakie
import PhasedChargeWatershed: load_data, build_viewer

isempty(ARGS) && error("usage: julia --project viewer.jl <data.json>")
pd = load_data(ARGS[1])

# the set of open windows,the program runs until the last one is closed
screens = Set{GLMakie.Screen}()

function add_view(init = nothing)
    fig = build_viewer(pd; on_add_view = add_view, init = init)
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
