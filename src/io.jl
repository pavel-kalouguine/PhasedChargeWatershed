
function load_data(file_path::String)::PhasedData
    data = parsefile(file_path)

    gens = [eval(Meta.parse(s)) for s in data["space_group"]]
    G = SpaceGroupQuotient(gens)
    @assert length(G) == length(gens) "error, elements do not form closed group"

    m = eval(Meta.parse(data["metric"]))
    D, N = size(m)
    md = SMatrix{D,N,Float64}(m)
    peaks=PhasedPeak{N}[]

    for r in data["reflections"]
        length(r["k"]) == N || throw(ArgumentError("reflection k has length $(length(r["k"])), expected $N"))
        k = SVector{N,Int}(r["k"])
        f_saved = Complex(r["ampl"][1], r["ampl"][2]) # The saved structure factor corresponding to `k`
        orbit=make_orbit(k, G)
        if orbit isa ExtinctOrbit
            throw(ArgumentError("The wavevector $k belongs to an extinct orbit."))
        end

        ϕ = first(ap.ϕ for ap in orbit.aps if ap.k == k) # Find the phase corresponding to `k` in the orbit.
        corr = exp(-2π * im * ϕ) # Phase factor needed to obtain the correct phase of the saved representative peak
        for ap in orbit.aps
            f = f_saved * corr * exp(2π * im * ap.ϕ)
            push!(peaks, PhasedPeak(ap.k, f))
            if orbit isa ComplexOrbit
                push!(peaks, PhasedPeak(-ap.k, conj(f))) # For orbits of complex type, add antipodes
            end
        end
    end

    PhasedData(G, md, peaks)
end