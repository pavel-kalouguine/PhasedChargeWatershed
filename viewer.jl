# Command-line launcher for the interactive density viewer.
#
#   julia --project viewer.jl path/to/data.json
#
# Loads a phased-data JSON file, opens the viewer window (GLMakie backend) and keeps
# it open until the user closes it.

using GLMakie
import PhasedChargeWatershed: load_data, build_viewer

isempty(ARGS) && error("usage: julia --project viewer.jl <data.json>")

pd  = load_data(ARGS[1])
fig = build_viewer(pd)

screen = display(fig)   # open the window
wait(screen)            # block until the user closes it
