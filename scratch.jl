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
    D, N = size(m)
    md = SMatrix{D,N,Float64}(m)
    dd = DiffractionData(G, md)
    sf = Vector{Complex{Float64}}(undef, length(data["reflections"]))
    for (i, r) in enumerate(data["reflections"])
        k = SVector{length(r["k"]),Int}(r["k"])
        orbit_length=add_peak!(dd, k, Float64(r["I"]))
        if orbit_length==0
            throw(ArgumentError("Reflection $k is already accounted for in the diffraction data"))
        end
        sf[i] = Complex(r["ampl"][1], r["ampl"][2])
    end
    # TODO: add a sanity check that the phases of structure factors of orbits of the real type are 
    # consistent with the symmetry of the space group. Attention: such structure factors are 
    # not necessarity real numbers!

    dd, sf
end

filepath = joinpath(@__DIR__, "data", "synthetic.json")
dd, sf = load_data(filepath)
