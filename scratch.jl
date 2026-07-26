using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density, create_watershed_grid, PhasedData, PhasedPeak, WatershedGrid, pre_watershed, WatershedResult, sample_pre_watershed_labels, quadratic_moment_on_neighbors



filepath = joinpath(@__DIR__, "data", "synthetic_pg.json")
phased_data = load_data(filepath)

result=pre_watershed(phased_data, density_factor=10.0)
sg=SamplingGrid{2,2}(SMatrix{2,2,Int}([1 0; 0 1]), SVector{2,Float64}(0.0, 0.0), (1024,1024))
sampled_labels=sample_pre_watershed_labels(result, sg)