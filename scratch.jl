using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density



filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)

grid=SamplingGrid{2,2}(SMatrix{2,2,Int}([1 0; 0 1]), SVector{2,Float64}(0.0, 0.0), (1024, 1024))
ρ=sample_density(phased_data.peaks, grid)
extrema(ρ)