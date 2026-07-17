


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
    ρ = zeros(Float64, grid.size...)
    sample_density!(ρ, peaks, grid)
    return ρ
end


function sample_density!(
    ρ::Array{Float64,M},
    peaks::AbstractVector{<:PhasedPeak{N}},
    grid::SamplingGrid{N,M}
)::Nothing where {N,M}

    @assert size(ρ) == grid.size "dimension mismatch between ρ and grid"

    h = div(grid.size[1], 2) + 1
    half_size = (h, grid.size[2:end]...)

    # 1. Allocate the complex accumulator locally
    fp = zeros(Complex{Float64}, half_size...)

    # 2. Compute the Fourier coefficients
    for p in peaks
        kp = alias(grid.direction * p.k, grid.size)
        if kp[1] <= h
            fp[kp...] += p.f * exp(1im * 2π * (p.k ⋅ grid.origin))
        end
    end

    # 3. Plan the IRFFT and execute directly into the user-provided ρ array
    plan = plan_irfft(fp, grid.size[1])
    mul!(ρ, plan, fp)

    return nothing
end