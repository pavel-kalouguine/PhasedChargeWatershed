using GLMakie, StaticArrays

include(joinpath(@__DIR__, "load_data.jl"))
include(joinpath(@__DIR__, "sampling.jl"))

filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)

grid = SamplingGrid(SA[2 0; 0 2], SA[0.5, 0.5], (1024, 1024))
ρ = sample_density(phased_data.peaks, grid)

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
heatmap!(ax, ρ)

fig