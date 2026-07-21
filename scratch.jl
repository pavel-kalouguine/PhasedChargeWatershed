using GLMakie, StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density


filepath = joinpath(@__DIR__, "data", "CdYb2.json")
phased_data = load_data(filepath)

grid = SamplingGrid(SA[0 1 1 1 1 -1; 1 0 0 0 0 0], SA[0., 0., 0., 0., 0., 0.], (843, 377))
ρ = sample_density(phased_data.peaks, grid)

fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
heatmap!(ax, ρ)

fig