


# TODO: consider using FFTViews.jl utility instead of this function
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
- performs an inverse Fourier transform using `irfft`  
- and returns the inverse real density values on the sampling grid.

Only the non-redundant half-spectrum along the first axis is stored, consistent
with `irfft` input layout.
"""
function sample_density(peaks::AbstractVector{<:PhasedPeak{N}}, grid::SamplingGrid{N,M})::Array{Float64,M} where {N,M}
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