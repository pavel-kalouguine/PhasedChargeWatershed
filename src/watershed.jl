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
        Q .+= abs2(f) * (k * k')
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
    sqrtQ=sqrt(Q) # Q does not change across attempts, so its square root is computed once
    best_candidate=nothing
    smallest_val=Inf
    for i in 1:n_attempts
        R = svd(randn(N, N)).U # Random orthogonal matrix
        B=R*B0 # Rotated basis
        L=round.(Int, scaling_factor*(B\sqrtQ))
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
    end
    best_candidate
end


"""
    pre_watershed(phased_data::PhasedData{N,D}; density_factor::Float64=1.0, n_attempts::Int=1000) where {N,D}

Perform the pre-watershed segmentation of the density sampled from phased peak data.

Constructs an optimal `WatershedGrid` via `create_watershed_grid`, samples the electron
density onto it, and then runs a watershed algorithm by processing grid sites in
decreasing order of density. Each site is assigned a basin label according to its
highest-density labeled neighbor; sites with no labeled neighbors seed new basins and
are recorded as summits. Saddle points between distinct basins are detected and stored.
Initially, every summit forms its own basin; basins are to be merged in a subsequent
post-processing step.

# Arguments
- `phased_data`: the phased peak data used to construct the watershed grid and sample
  the density.
- `density_factor`: passed to `create_watershed_grid`; must be `>= 1.0` and controls
  the density of grid sites relative to a heuristic minimum.
- `n_attempts`: number of random candidate lattices tried by `create_watershed_grid`.

Returns a `WatershedResult` containing the phased data itself, the sampled density,
basin labels, summit indices, saddle points, and basin assignment array.

Throws an `ArgumentError` if no valid `WatershedGrid` could be constructed within
`n_attempts` attempts.
"""
function pre_watershed(phased_data::PhasedData{N,D}; density_factor::Float64=1.0, n_attempts::Int=1000) where {N,D}
    n_attempts >= 1 || throw(ArgumentError("n_attempts must be >= 1"))
    wg=create_watershed_grid(phased_data; density_factor=density_factor, n_attempts=n_attempts)
    if wg === nothing
        throw(ArgumentError("Failed to create a valid WatershedGrid after $n_attempts attempts."))
    end
    # Pre-allocate the watershed result
    result=WatershedResult(phased_data, wg)
    sample_density!(result.ρ, phased_data.peaks, wg.grid)
    
    # Sort the indices of grid sites by decreasing density
    sorted_indices=sortperm(result.ρ, rev=true)
    # Pre-allocate arrays to store the labels and densities of the neighbors of each site
    labeled_neighbors=zeros(Int, length(wg.neighbors)) 
    density_of_labeled_neighbors=zeros(Float64, length(wg.neighbors)) 
    shifts_of_labeled_neighbors=zeros(Int, length(wg.neighbors))
    d=wg.grid.size[1]
    for i in sorted_indices
        # Look for the labeled neighbors
        more_than_one_neighboring_label=false # Are we at the watershed line?
        num_labeled_neighbors=0
        last_found_label=0
        highest_density_of_labeled_neighbors=-Inf
        label_with_highest_density=0
        for shift in wg.neighbors
            j=mod1(i+shift, d)
            label=result.labels[j]
            if label != 0
                if result.ρ[j] > highest_density_of_labeled_neighbors
                    highest_density_of_labeled_neighbors=result.ρ[j] # Keep track of the highest density among the labeled neighbors
                    label_with_highest_density=label
                end
                num_labeled_neighbors+=1
                labeled_neighbors[num_labeled_neighbors]=label
                density_of_labeled_neighbors[num_labeled_neighbors]=result.ρ[j]
                shifts_of_labeled_neighbors[num_labeled_neighbors]=shift
                if last_found_label != 0 && label != last_found_label
                    more_than_one_neighboring_label=true
                end
                last_found_label=label
            end
        end
        # If there are no labeled neighbors, this is a new basin
        if num_labeled_neighbors == 0
            push!(result.summits, i)
            result.labels[i]=length(result.summits) # Assign a new label to the current site        
        else 
            # Assign to the site the label of the neighbor with the highest density
            result.labels[i]=label_with_highest_density
        end
        if more_than_one_neighboring_label
            # If there are multiple neighboring labels, there is one or more saddle points. 
            for n=1:num_labeled_neighbors
                labeled_neighbors[n] == result.labels[i] && continue # Skip neighbors sharing the site's own label; not a saddle
                sp=SaddlePoint(
                    (result.labels[i], labeled_neighbors[n]),
                    (i, mod1(i+shifts_of_labeled_neighbors[n], d)),
                    (result.ρ[i], density_of_labeled_neighbors[n])
                ) 
                if !haskey(result.saddles, sp.labels)
                    result.saddles[sp.labels]=sp
                end   
            end
        end 
    end
    resize!(result.basins, length(result.summits))
    result.basins.= 1:length(result.summits) # Every summit has its own basin, initially. The basins will be merged in the post-processing stage.
    result
end



"""
    basin_root(parent::Vector{Int}, i::Int)

Index of the basin the pre-basin `i` was fused into, following `parent` up to a fixed point.
"""
function basin_root(parent::Vector{Int}, i::Int)
    while parent[i] != i
        i = parent[i]
    end
    return i
end


"""
    fuse_basins!(parent::Vector{Int}, i::Int, j::Int)

Fuse the basins of the pre-basins `i` and `j`, keeping the smaller of the two indices.
"""
function fuse_basins!(parent::Vector{Int}, i::Int, j::Int)
    a, b = minmax(basin_root(parent, i), basin_root(parent, j))
    parent[b] = a
    return nothing
end


"""
    update_basins!(result::WatershedResult, summit_percent::Real, saddle_ratio::Real, reference::Function)

Recompute `result.basins` in place.

Pre-basins are fused when the saddle point between them rises above `saddle_ratio` of their
summits, compared as `reference` says (`min`, `max`, or the average); the fused group takes
the smallest of its indices, which is also its highest summit. A group whose summit stays
below `summit_percent` of the range the summits span is left unlabeled, its elements set
to 0. The range is taken between the lowest and the highest summit rather than from zero:
the density is only known up to a constant, so zero is not a physical level.

`saddle_ratio = 1` together with `summit_percent = 0` reproduces the pre-watershed
labelling exactly: a saddle point never rises above the summits it separates, and the
interpolation of the cutoff is exact at the ends of the range.
"""
function update_basins!(result::WatershedResult, summit_percent::Real, saddle_ratio::Real,
                        reference::Function)
    summit_ρ = result.ρ[result.summits]
    parent = collect(eachindex(summit_ρ))
    for sp in values(result.saddles)
        a, b = sp.labels
        if minimum(sp.values) > saddle_ratio * reference(summit_ρ[a], summit_ρ[b])
            fuse_basins!(parent, a, b)
        end
    end
    lowest, highest = extrema(summit_ρ)
    t = summit_percent / 100
    cutoff = (1 - t) * lowest + t * highest
    for i in eachindex(summit_ρ)
        r = basin_root(parent, i)
        result.basins[i] = summit_ρ[r] < cutoff ? 0 : r
    end
    return result
end


"""
    sample_pre_watershed_labels(result::WatershedResult{N}, grid::SamplingGrid{N,M}) where {N,M}

Sample the pre-watershed basin labels onto the sites of a given sampling grid.

For each site in `grid`, the corresponding position in the `N`-dimensional unit cell is
computed, and a corresponding site index in the cyclic watershed grid is estimated by
rounding in the watershed lattice basis (see the TODO in the implementation).

Returns an `Array{Int,M}` of size `grid.size`, where each element is the integer basin
label of the watershed grid site nearest to the corresponding sampling point of `grid`.
"""
function sample_pre_watershed_labels(result::WatershedResult{N}, grid::SamplingGrid{N,M})::Array{Int,M} where {N,M}
    d=result.wg.grid.size[1]
    sampled_labels=zeros(Int, grid.size...)
    # Loop over the sites of the sampling grid
    for ind in CartesianIndices(grid.size)
        # Position of the site in the `N`-dimensional unit cell
        x = mod.(grid.direction' * SVector((ind.I.-1)./grid.size) + grid.origin, 1) 
        v=round.(Int, result.wg.L * x) # Nearest site position in the `N`-dimensional basis of the watershed grid
        #TODO: Find the actual nearest site
        # Index of the nearest site in the cyclic watershed grid. The shift given by the
        # basis indices is counted from site 1, which sits at the origin.
        i = mod(v'*result.wg.basis_indices, d) + 1
        sampled_labels[ind] = result.labels[i]
    end
    sampled_labels
end