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
struct SamplingGrid{N,M}
    direction::SMatrix{M,N,Int}
    origin::SVector{N,Float64}
    size::NTuple{M,Int}
end

# Copy-pasted from ChargeFlipPhaser, probably should export it properly
alias(k::SVector{M,Int}, size::NTuple{M,Int}) where M =
    tuple((mod.(k, SVector{M,Int}(size...)) .+ 1)...)


"""
    sample_density(peaks, grid) -> Array{Float64,M}

Compute the real-space density on the affine sampling lattice described by `grid`
from phased Fourier peaks `peaks`.

The function:
- projects each peak wavevector onto the sampling lattice using `grid.direction`,
- applies periodic indexing with `alias`,
- accumulates complex Fourier coefficients including the origin phase shift
  `exp(1im * 2π * p.k ⋅ grid.origin)`,
- and returns the inverse real FFT (`irfft`) over the first dimension.

Only the non-redundant half-spectrum along the first axis is stored, consistent
with `irfft` input layout.
"""
function sample_density(peaks::Vector{PhasedPeak{N}}, grid::SamplingGrid{N,M})::Array{Float64,M} where {N,M}
    h=div(grid.size[1], 2) + 1 # Half of the size of the first dimension, for irfft
    # Compute the Fourier coefficients
    half_size = (h, grid.size[2:end]...)
    fp=zeros(Complex{Float64}, half_size...) # Accumulator array
    for p in peaks
        kp=alias(grid.direction * p.k, grid.size) # Projected wavevector
        if kp[1] <= h  # Only accumulate the first half of the first dimension
            fp[kp...] += p.f*exp(1im * 2π * p.k ⋅ grid.origin) # Take the origin into account
        end
    end
    irfft(fp, grid.size[1]) # Inverse Fourier transform to get the density
end