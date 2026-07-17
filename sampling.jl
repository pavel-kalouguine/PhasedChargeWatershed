using StaticArrays
using FFTW
using LinearAlgebra

"""
        SamplingGrid{N, D}

Parameters of a `M`-dimensional affine lattice commensurate to an `N`-dimensional
periodic lattice.

- `direction`: rows are the `M` vectors of the `N`-dimensional lattice common
    to the affine sublattice.
- `size`: integer factors such that each row in `direction` is a multiple of a
    basis vector of the affine lattice.
- `origin`: translation of the affine lattice origin with respect to the
    origin of the `N`-dimensional lattice.

This structure can be used both for the watershed algorithm and for the visualization 
(with `M=2`)
"""
struct SamplingGrid{N, M}
    direction::SMatrix{M,N, Int}
    origin::SVector{N,Float64}
    size::NTuple{M, Int}
end

# Copy-pasted from ChargeFlipPhaser, probably should export it properly
alias(k::SVector{M,Int}, size::NTuple{M,Int}) where M =
    tuple((mod.(k, SVector{M,Int}(size...)) .+ 1)...)
 
# TODO: use irfft instead of complex ifft for speed
function sample_density(peaks::Vector{PhasedPeak{N}}, grid::SamplingGrid{N, M})::Array{Float64, M} where {N, M}
    # Compute the Fourier coefficients
    fp=zeros(Complex{Float64}, grid.size...) # Accumulator array
    for p in peaks
        kp=alias(grid.direction * p.k, grid.size) # Projected wavevector
        fp[kp...]+= p.f*exp(1im * 2π * p.k ⋅ grid.origin) # Take the origin into account
    end
    real.(ifft(fp))
end