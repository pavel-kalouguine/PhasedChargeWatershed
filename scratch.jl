using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density, create_watershed_grid, PhasedData, PhasedPeak, WatershedGrid



filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)

wg=create_watershed_grid(phased_data)