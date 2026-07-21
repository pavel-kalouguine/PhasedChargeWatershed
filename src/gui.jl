# Interactive viewer for a phased density.
#
# The window shows a 2D section of the density on top, and a single row of controls
# underneath: the sampling grid's
#   - `direction` : a 2×N integer matrix (which two lattice directions span the cut),
#   - `origin`    : an N-vector, each component in [0, 1] (where the cut sits),
#   - `size`      : the (nx, ny) resolution of the output image.
# N is the number of indices in each wave vector and comes from the data, so the
# controls are built to fit whatever N the file has.
#
# The code uses only the Makie API (not a concrete backend), so the launcher picks
# the backend (GLMakie for an interactive window).

"""
    build_viewer(pd::PhasedData) -> Figure

Build an interactive window that shows a 2D section of the density stored in `pd`
and lets the user change the sampling grid (`direction`, `origin`, `size`) live.

On open it shows a default section (the first two lattice directions, origin at 0,
1024×1024); the user can then edit any parameter to explore other sections.
"""
function build_viewer(pd::PhasedData{N}) where {N}
    fig = Figure(size = (950, 950))

    #density image
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = "density section")
    hidedecorations!(ax)

    #controls
    controls = GridLayout(fig[2, 1], tellheight = true)

    #direction
    dircol = GridLayout(controls[1, 1])
    Label(dircol[1, 1:N], "direction (2×$N)")
    dir_boxes = [Textbox(dircol[1+i, j]; width = 46, validator = Int,
                         stored_string = string(i == j ? 1 : 0))
                 for i in 1:2, j in 1:N]

    # origin
    origcol = GridLayout(controls[1, 2])
    Label(origcol[1, 1:N], "origin (0…1)")
    orig_sliders = [Slider(origcol[2, j]; range = 0:0.01:1, startvalue = 0.0, width = 90) for j in 1:N]
    orig_boxes   = [Textbox(origcol[3, j]; width = 60, validator = Float64, stored_string = "0.0") for j in 1:N]
    for j in 1:N
        # slider moved = show the value in its box
        on(orig_sliders[j].value) do v
            orig_boxes[j].displayed_string[] = string(round(v; digits = 2))
        end
        
        on(orig_boxes[j].stored_string) do s
            v = tryparse(Float64, something(s, ""))
            v === nothing || set_close_to!(orig_sliders[j], v)
        end
    end

    # size
    sizecol = GridLayout(controls[1, 3])
    Label(sizecol[1, 1:2], "size")
    size_boxes = [Textbox(sizecol[2, j]; width = 70, validator = Int, stored_string = "1024") for j in 1:2]

    
    read_int(tb, default) = something(tryparse(Int, something(tb.stored_string[], "")), default)

    #Assemble the current sampling grid from all the controls.
    function current_grid()
        d = SMatrix{2,N,Int}([read_int(dir_boxes[i, j], i == j ? 1 : 0) for i in 1:2, j in 1:N])
        o = SVector{N,Float64}([sl.value[] for sl in orig_sliders])
        s = (max(2, read_int(size_boxes[1], 1024)), max(2, read_int(size_boxes[2], 1024)))
        SamplingGrid(d, o, s)
    end

    grid = Observable(current_grid())
    for tb in dir_boxes;    on(_ -> (grid[] = current_grid()), tb.stored_string); end
    for sl in orig_sliders; on(_ -> (grid[] = current_grid()), sl.value);         end
    for tb in size_boxes;   on(_ -> (grid[] = current_grid()), tb.stored_string); end

    density = lift(g -> sample_density(pd.peaks, g), grid)

    heatmap!(ax, density; colormap = :jet,
             colorrange = lift(ρ -> (minimum(ρ), maximum(ρ) + eps()), density))

    fig
end
