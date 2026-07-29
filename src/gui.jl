# Interactive viewer for a phased density.
#
# The window shows a 2D section of the density on top (with a colorbar to its right),
# and a single row of controls underneath: the sampling grid's
#   - `direction` : a 2×N integer matrix (which two lattice directions span the cut),
#   - `origin`    : an N-vector, each component in [0, 1] (where the cut sits),
#   - `size`      : the (nx, ny) resolution of the output image.
# N comes from the data, so the controls are built to fit whatever N the file has.
# The section is (re)computed when the "Apply" button is pressed, so heavy data are
# recomputed only on demand. The image aspect follows the requested size.
#
# When a watershed result is given, its basin boundaries are drawn over the density. The
# postprocessing of the basins is driven from a window of its own and acts on every view at
# once. It does not resample the section: it only rewrites which pre-basin belongs to which
# basin, so the boundaries follow the sliders as they are dragged.
#
# The code uses only the Makie API; the launcher picks the backend (GLMakie).

"""
    global_density_limits(pd::PhasedData) -> Tuple{Float64,Float64}

Estimate the global minimum and maximum of the density over the whole unit cell, for
use as a common colour range in all section heatmaps. The density is sampled once on a
minimal watershed grid (default `density_factor`), which covers the cell densely and
respects the data anisotropy; the extrema are slightly underestimated (finite grid), so
a few points near the true min/max may saturate. Returns `nothing` if no grid was built.
"""
function global_density_limits(pd::PhasedData)
    wg = create_watershed_grid(pd)
    wg === nothing && return nothing
    ρ = sample_density(pd.peaks, wg.grid)
    lo, hi = extrema(ρ)
    return (lo, hi > lo ? hi : hi + eps(hi))
end

"""
    global_density_limits(result::WatershedResult) -> Tuple{Float64,Float64}

Same limits taken from an existing `WatershedResult`, whose density has already been
sampled over the whole cell. No grid is built, so this is immediate.
"""
function global_density_limits(result::WatershedResult)
    lo, hi = extrema(result.ρ)
    return (lo, hi > lo ? hi : hi + eps(hi))
end

"""
    basin_boundaries(labels::Array{Int,2}, width::Int = 2) -> Matrix{Float32}

Mark the sites lying on a boundary between two different basins, drawing a line `width`
sites across. Whenever two adjacent sites carry different labels, the line is laid down
symmetrically about the pair, so it straddles the boundary instead of sitting on one side
of it. The cut is periodic by design, so the neighbours are taken modulo the size of the
image. Sites inside a basin get `NaN`, so that the result can be drawn as a heatmap whose
`nan_color` is transparent, leaving the density below it visible.
"""
function basin_boundaries(labels::Array{Int,2}, width::Int = 2)
    nx, ny = size(labels)
    mask = fill(NaN32, nx, ny)
    offsets = (1 - width ÷ 2):(width - width ÷ 2)
    for i in 1:nx, j in 1:ny
        l = labels[i, j]
        if labels[mod1(i + 1, nx), j] != l
            for d in offsets
                mask[mod1(i + d, nx), j] = 1.0f0
            end
        end
        if labels[i, mod1(j + 1, ny)] != l
            for d in offsets
                mask[i, mod1(j + d, ny)] = 1.0f0
            end
        end
    end
    return mask
end

"""
    boundary_width(image_size) -> Int

How many sites wide the boundary lines should be drawn. A fixed number of sites is not a
fixed thickness on screen: the section is scaled to the window, so a two-site line that
reads well on a 512×512 cut falls below one screen pixel on a 2048×2048 one. Growing the
width with the section keeps the lines about equally thick whatever the resolution.
"""
boundary_width(image_size) = max(2, round(Int, maximum(image_size) / 256))

"""
    build_viewer(pd::PhasedData; on_add_view = _ -> nothing, init = nothing, colorrange = nothing, result = nothing) -> Figure

Build an interactive window that shows a 2D section of the density stored in `pd`.
The sampling grid (`direction`, `origin`, `size`) is edited in the controls and
applied with the "Apply" button. "Add a view" calls `on_add_view(grid)` with the
current sampling grid, so the caller can open another window; `init` is a `SamplingGrid`
the controls start from (used to clone a view). On open it shows `init`, or the default
section (first two axes, origin 0, 1024×1024) when `init` is not given. `colorrange`, if
given, fixes the heatmap colour range for every section (a shared/global scale);
otherwise each section is scaled to its own extrema. `result`, if given, is an
`Observable` holding a `WatershedResult`, whose basin boundaries are drawn as white lines
on top of the density. Several windows can share the same observable, and then they all
follow the postprocessing driven from `build_basin_controls`.
"""
function build_viewer(pd::PhasedData{N}; on_add_view = _ -> nothing, init = nothing, colorrange = nothing, result = nothing) where {N}
    init_grid = init === nothing ?
        SamplingGrid(SMatrix{2,N,Int}([i == j ? 1 : 0 for i in 1:2, j in 1:N]),
                     zero(SVector{N,Float64}), (1024, 1024)) : init

    fig = Figure(size = (950, 950))

    #density image with a colorbar glued to its right
    top = GridLayout(fig[1, 1])
    got_global = colorrange !== nothing
    scale_note = got_global ? "global colour scale" : "per-section colour scale (global density limits unavailable)"
    ax = Axis(top[1, 1], aspect = DataAspect(), title = "density section",
              subtitle = scale_note, subtitlecolor = got_global ? :black : :red)
    hidedecorations!(ax)

    #controls
    controls = GridLayout(fig[2, 1], tellheight = true)

    #direction
    dircol = GridLayout(controls[1, 1])
    Label(dircol[1, 1:N], "direction (2×$N)")
    step_box!(tb, d) = (tb.displayed_string[] =
        string(something(tryparse(Int, something(tb.displayed_string[], "")), 0) + d))
    dir_boxes = Matrix{Textbox}(undef, 2, N)
    for i in 1:2, j in 1:N
        cell = GridLayout(dircol[1+i, j])
        tb = Textbox(cell[1:2, 1]; width = 38, validator = Int, stored_string = string(init_grid.direction[i, j]))
        up = Button(cell[1, 2]; label = "▲", width = 18, fontsize = 8)
        dn = Button(cell[2, 2]; label = "▼", width = 18, fontsize = 8)
        on(_ -> step_box!(tb, 1), up.clicks)
        on(_ -> step_box!(tb, -1), dn.clicks)
        dir_boxes[i, j] = tb
    end

    #origin
    origcol = GridLayout(controls[1, 2])
    Label(origcol[1, 1:N], "origin (0…1)")
    orig_sliders = [Slider(origcol[2, j]; range = 0:0.01:1, startvalue = init_grid.origin[j], width = 90) for j in 1:N]
    orig_boxes   = [Textbox(origcol[3, j]; width = 60, validator = Float64, stored_string = string(init_grid.origin[j])) for j in 1:N]
    for j in 1:N
        on(orig_sliders[j].value) do v
            orig_boxes[j].displayed_string[] = string(round(v; digits = 2))
        end
        on(orig_boxes[j].stored_string) do s
            v = tryparse(Float64, something(s, ""))
            v === nothing || set_close_to!(orig_sliders[j], v)
        end
    end

    #size
    sizecol = GridLayout(controls[1, 3])
    Label(sizecol[1, 1:2], "size")
    size_boxes = [Textbox(sizecol[2, j]; width = 70, validator = Int, stored_string = string(init_grid.size[j])) for j in 1:2]

    #apply / add a cloned view
    apply   = Button(controls[1, 4], label = "Apply")
    addview = Button(controls[1, 5], label = "Add a view")

    function read_int(tb, default)
        s = something(tb.displayed_string[], something(tb.stored_string[], ""))
        something(tryparse(Int, s), default)
    end

    # Assemble the sampling grid from all the controls.
    function current_grid()
        d = SMatrix{2,N,Int}([read_int(dir_boxes[i, j], init_grid.direction[i, j]) for i in 1:2, j in 1:N])
        o = SVector{N,Float64}([sl.value[] for sl in orig_sliders])
        s = (max(2, read_int(size_boxes[1], init_grid.size[1])), max(2, read_int(size_boxes[2], init_grid.size[2])))
        SamplingGrid(d, o, s)
    end

    # Recompute only when Apply is pressed, and only if the controls really describe
    # another cut: pressing Apply again with the same settings would otherwise redo the
    # whole sampling and throw the identical result away.
    grid = Observable(current_grid())
    on(apply.clicks) do _
        new_grid = current_grid()
        new_grid == grid[] || (grid[] = new_grid)
    end
    on(_ -> on_add_view(grid[]), addview.clicks)

    density = lift(g -> sample_density(pd.peaks, g), grid)

    cr = colorrange === nothing ? lift(ρ -> (minimum(ρ), maximum(ρ) + eps()), density) : colorrange
    hm = heatmap!(ax, density; colormap = :jet, colorrange = cr)

    Colorbar(top[1, 2], hm, label = "density")
    colsize!(top, 1, Aspect(1, init_grid.size[1] / init_grid.size[2]))

    if result !== nothing
        prelabels = lift(g -> sample_pre_watershed_labels(result[], g), grid)
        boundaries = lift(prelabels, result) do pl, r
            basin_boundaries(map(l -> r.basins[l], pl), boundary_width(size(pl)))
        end
        heatmap!(ax, boundaries; colormap = [:white, :white], colorrange = (0, 1),
                 nan_color = :transparent)
    end

    on(_ -> reset_limits!(ax), density)

    resize_to_layout!(fig)

    fig
end

"""
    build_basin_controls(result::Observable{<:WatershedResult}) -> Figure

Build the window postprocessing the basins of `result`.

The postprocessing applies to the `N`-dimensional basins, of which the viewers only show
2D cuts, so it is driven from a window of its own and acts on every open view at once.
"""
function build_basin_controls(result)
    fig = Figure(size = (600, 170))

    summit = Slider(fig[1, 2]; range = 0:0.001:1, startvalue = 0, width = 260)
    Label(fig[1, 1], lift(summit.value) do v
              "unlabel basins whose summit is below $(round(100v; digits = 1))% of the summit range"
          end, halign = :right)
    ratio = Slider(fig[2, 2]; range = 0:0.005:1, startvalue = 1, width = 260)
    Label(fig[2, 1], lift(ratio.value) do v
              "fuse basins whose saddle is above $(round(100v; digits = 1))% of their summits"
          end, halign = :right)
    references = [("average", ((a, b) -> (a + b) / 2)), ("minimum", min), ("maximum", max)]
    reference = Menu(fig[3, 2]; options = references, default = "average", width = 260)
    Label(fig[3, 1], "summits compared as their", halign = :right)

    for signal in (summit.value, ratio.value, reference.selection)
        on(signal) do _
            update_basins!(result[], summit.value[], ratio.value[], reference.selection[])
            notify(result)
        end
    end

    resize_to_layout!(fig)

    fig
end
