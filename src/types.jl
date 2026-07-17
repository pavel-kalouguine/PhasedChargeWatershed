"""
    PhasedPeak{N}

Represents a phased Bragg peak.

The data stores the wave vector and the complex structure factor.
"""
struct PhasedPeak{N}
    k::SVector{N,Int}
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
struct PhasedData{N,D}
    G::SpaceGroupQuotient{N}
    md::SMatrix{D,N,Float64}
    peaks::Vector{PhasedPeak{N}}
end

"""
        SamplingGrid{N, M}

Parameters of a `M`-dimensional affine lattice commensurate to an `N`-dimensional
periodic lattice.

- `direction`: rows are the `M` vectors of the `N`-dimensional lattice common
    to the affine sublattice.
- `size`: number of sampling points along each of the `M` affine lattice directions.
     This corresponds to the output array shape.
- `origin`: translation of the affine lattice origin with respect to the
    origin of the `N`-dimensional lattice.

This structure can be used both for the watershed algorithm and for the visualization 
(with `M=2`)
"""
struct SamplingGrid{N,M}
    direction::SMatrix{M,N,Int}
    origin::SVector{N,Float64}
    size::NTuple{M,Int}
end