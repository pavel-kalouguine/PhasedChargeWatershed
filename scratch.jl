using SpaceGroups
using ChargeFlipPhaser
using StaticArrays
using JSON


"""
    PhasedPeak{N}

Represents a phased Bragg peak.

The data store the wave vector and the complex structure factor.
"""
struct PhasedPeak{N}
    k::SVector{N, Int}
    f::Complex{Float64}
end



"""
        PhasedData{N, D}

Represents the result of the phasing.

Fields:
- `G` is the space group quotient.
- `md` is the metric data of the lattice; its columns are the basis vectors of the
    (quasi) lattice in the `D`-dimensional physical space.
- `peaks` is the collection of all phased peaks, including those related by the
    symmetry transformations.
"""
struct PhasedData{N, D}
    G::SpaceGroupQuotient{N}
    md::SMatrix{D,N,Float64}
    peaks::Vector{PhasedPeak{N}}
end


function load_data(file_path::String)::PhasedData
    data = JSON.parsefile(file_path)

    gens = [eval(Meta.parse(s)) for s in data["space_group"]]
    G = SpaceGroupQuotient(gens)
    @assert length(G) == length(gens) "error, elements do not form closed group"

    m = eval(Meta.parse(data["metric"]))
    D, N = size(m)
    md = SMatrix{D,N,Float64}(m)
    peaks=PhasedPeak{N}[]

    
    for (i, r) in enumerate(data["reflections"])
        k = SVector{length(r["k"]),Int}(r["k"]) # The wave vector (one per orbit)        
        f = Complex(r["ampl"][1], r["ampl"][2]) # The structure factor corresponding to `k`
        orbit=make_orbit(k, G)
        if orbit isa ExtinctOrbit
            throw(ArgumentError("The wavevector $k belongs to an extinct orbit."))
        end
        # TODO: Create a `PhasedPeak` for each element of the orbit (adding antipodes in the case of `RealOrbit`)
        # Carefully keep track of phase factors, the saved representative should enter with the factor as saved,
        # while the phases of the other elements of the orbiit should take into account the differences of the 
        # phases of the corresponding `AffinePhase` objects. In the case of `ComplexOrbit`, the structure factor 
        # of antipodal peaks should always be conjugates. 
    end

    PhasedData(G, md, peaks)
end

filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)
