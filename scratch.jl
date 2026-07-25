using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density, create_watershed_grid, PhasedData, PhasedPeak, WatershedGrid, pre_watershed, WatershedResult



filepath = joinpath(@__DIR__, "data", "synthetic_pg.json")
phased_data = load_data(filepath)

wg=create_watershed_grid(phased_data)
result=pre_watershed(phased_data)