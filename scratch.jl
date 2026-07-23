using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, create_watershed_grid



filepath = joinpath(@__DIR__, "data", "CdYb.json")
phased_data = load_data(filepath)



create_watershed_grid(phased_data)