include(joinpath(@__DIR__, "load_data.jl"))

filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)
