"""
    quadratic_moment(peaks::Vector{PhasedPeak{N}}) where N

Compute the covariance matrix (quadratic moment) of a collection of Bragg peaks.

For each peak with reciprocal lattice vector `k` and complex structure factor `f`,
the contribution `abs2(f) * (k * k')` is accumulated, and the result is normalized
by the total intensity `sum(abs2(f))` over all peaks. This yields the intensity-weighted
second moment `⟨k k'⟩` of the peak positions, i.e. the covariance matrix of the
distribution of reciprocal lattice vectors weighted by their intensities.

Returns a `Symmetric{Float64}` `N x N` matrix.
"""
function quadratic_moment(peaks::Vector{PhasedPeak{N}}) where N
    # Compute the quadratic moment of the phased peaks
    Q = zeros(Float64, N, N)
    for peak in peaks
        k = peak.k
        f = peak.f
        Q += abs2(f) * (k * k')
    end
    total_intensity = sum(abs2.(getfield.(peaks, :f)))
    if total_intensity == 0.0
        throw(ArgumentError("Total intensity of peaks is zero; cannot compute quadratic moment."))
    end
    return Symmetric(Q/total_intensity)
end


"""
    quadratic_moment_on_neighbors(wg::WatershedGrid{N}) where N

Compute, for each neighbor of the origin in the watershed grid `wg`, the value of the 
quadratic form `wg.Q` evaluated at the corresponding fractional lattice site.

The watershed grid forms a cyclic lattice, which can be indexed by a single integer.
The position of the neighbor indexed by `shift` is a rational vector `x` computed as
`x = direction * shift / d`, wrapped into the range `[-1/2, 1/2]` by
subtracting the nearest integer (`round.(x)`), where `direction` is the
grid's direction matrix and `d` is the grid size along its first axis (the actual 
computation is performed with floating point numbers to avoid integer overflow).
The returned value for that neighbor is `x' * wg.Q * x`.

Returns a `Vector{Float64}` with one entry per neighbor in `wg.neighbors`.
"""
function quadratic_moment_on_neighbors(wg::WatershedGrid{N}) where N
    if size(wg.Q) != (N, N)
        throw(ArgumentError("Quadratic moment must be an $N x $N matrix"))
    end
    grid = wg.grid
    direction = grid.direction
    dir=SVector(direction')
    d=grid.size[1]
    neighbors = wg.neighbors
    vals=Vector{Float64}(undef, length(neighbors))
    for n in eachindex(neighbors)
        shift=neighbors[n]
        x=dir*shift/d
        x=x-round.(x)
        vals[n]=x' * wg.Q * x
    end
    return vals
end


"""
    create_watershed_grid(phased_data::PhasedData{N,D}; density_factor::Float64=1.0, n_attempts::Int=1000) where {N,D}

Construct a `WatershedGrid` for the given phased data, optimized for the
anisotropy of the underlying reciprocal-space intensity distribution.

The intensity-weighted covariance (quadratic moment) `Q` of the Bragg peaks
is used to characterize the anisotropy of the data. The neighbor topology of
the watershed grid is fixed to match the neighbor graph of a known dense
sphere packing in dimension `N` (i.e. the graph whose edges connect sites at
the shortest packing distance). Within this fixed topology, many random
cyclic sublattices realizing that neighbor structure are generated (via
random rotations of the packing basis followed by Smith normal form
reduction to a cyclic grid), and for each candidate grid the values of the
quadratic form `Q` on the vectors linking neighboring sites are computed
using `quadratic_moment_on_neighbors`. The candidate minimizing the maximum
such value (i.e. the one whose neighbor edges are as short as possible with
respect to the data's covariance metric) is kept.

# Arguments
- `phased_data`: the phased peak data, from which the covariance quadratic
  form `Q` is computed via `quadratic_moment`.
- `density_factor`: must be `>= 1.0`; controls the overall scaling of the
  lattice, and hence the density of grid sites, relative to a heuristic
  minimal scaling factor.
- `n_attempts`: number of random candidate lattices to try; the best one
  found (smallest maximal quadratic form value on neighbor edges) is
  returned.

Returns the best `WatershedGrid` found, or `nothing` if no valid candidate
was found within `n_attempts` tries.
"""
function create_watershed_grid(phased_data::PhasedData{N,D}; density_factor::Float64=1.0, n_attempts::Int=1000) where {N,D}
    if density_factor <1.0
        throw(ArgumentError("density_factor must be >= 1.0"))
    end
    Q = quadratic_moment(phased_data.peaks)
    B0=basis_of_dense_packing(N)
    nbs = nearest_neighbors(B0)
    # TODO: find the reasonable value of the minimal_scaling_factor and document the implications of this choice. 
    minimal_scaling_factor=20.0 # A heuristic scaling factor to make the grid dense enough to capture the density variations
    scaling_factor=density_factor^(1/N)*minimal_scaling_factor
    i=0
    best_candidate=nothing
    smallest_val=Inf
    while true
        i+=1
        R = svd(randn(N, N)).U # Random orthogonal matrix
        B=R*B0 # Rotated basis
        L=round.(Int, scaling_factor*(B\sqrt(Q)))
        s=nothing
        try
            # TODO: consider using BigInt to avoid integer overflow, but this will be slower.
            # The current implementation may fail for large scaling_factor.
            s = snf(L)
        catch err
            err isa InterruptException && rethrow()
            continue # SNF failed, try again
        end
        divisors=diag(s)
        if any(divisors[1:end-1] .!= 1)
            continue # We are looking for a cyclic lattice, try again
        end
        d=divisors[end]
        dir=s.V[:, end]
        direction = SMatrix{1, N, Int}(dir')
        grid=SamplingGrid{N,1}(direction, zero(SVector{N,Float64}), (d,))
        basis_indices=SVector{N,Int}(mod.(s.U[end,:], d))
        neighbors=Int[]
        for u in nbs
            shift=basis_indices'*u
            push!(neighbors, mod(shift, d))
        end
        candidate=WatershedGrid(grid, neighbors, basis_indices, SMatrix{N,N,Int}(L), Q)
        vals=quadratic_moment_on_neighbors(candidate)
        max_val=maximum(vals)
        if max_val < smallest_val
            smallest_val=max_val
            best_candidate=candidate
        end
        if i>=n_attempts
            break
        end
        
    end
    best_candidate
end