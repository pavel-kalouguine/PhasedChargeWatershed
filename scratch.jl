using SpaceGroups
using ChargeFlipPhaser
using StaticArrays
using JSON

function load_data(file_path::String)
    data = JSON.parsefile(file_path)

    gens = [eval(Meta.parse(s)) for s in data["space_group"]]
    G = SpaceGroupQuotient(gens)
    @assert length(G) == length(gens) "error, elements do not form closed group"

    m = eval(Meta.parse(data["metric"]))
    md = SMatrix{size(m, 1),size(m, 2),Float64}(m)
    dd = DiffractionData(G, md)
    sf = Vector{Complex{Float64}}(undef, length(data["reflections"]))
    for (i, r) in enumerate(data["reflections"])
        k = SVector{length(r["k"]),Int}(r["k"])
        add_peak!(dd, k, Float64(r["I"]))
        sf[i] = Complex(r["ampl"][1], r["ampl"][2])
    end

    dd, sf
end

filepath = joinpath(@__DIR__, "data", "CdYb.json")
dd, sf = load_data(filepath)
