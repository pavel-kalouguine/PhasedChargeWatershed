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


"""
        WatershedGrid{N}

Parameters of a cyclic sampling grid, i.e. a `SamplingGrid{N,M}` with `M=1`,
together with the neighborhood relation used by the watershed algorithm.

In cyclic grids the sites are indexed by a single integer. The field
`neighbors` is a collection of such indices corresponding to the sites
declared as neighbors of the origin. In general, two sites of the grid are
neighbors if the difference of their indices modulo the size of the grid
belongs to the collection `neighbors`. This relation of neighborhood is used
by the watershed algorithm when comparing the values of the density on
neighboring sites.

Fields:
- `grid`: the underlying `SamplingGrid{N,1}`.
- `neighbors`: the indices, modulo the grid size, defining the neighborhood
    relation between sites.
- `basis_indices`: the indices of the basis vectors inherited from the dense
    sphere packing in the `N`-dimensional space.
- `L`: the integer matrix transforming the sampling grid to the lattice of
    integer points in the `N`-dimensional space. In particular, the columns
    of `inv(L)` correspond to the grid positions indexed by `basis_indices`.
- `Q`: the covariance matrix (quadratic moment) of the Bragg peaks. The
    sampling grid is constructed so that the values of the quadratic form
    defined by `Q` on the neighbors of the origin are close to each other.
"""
struct WatershedGrid{N}
    grid::SamplingGrid{N,1}    
    neighbors::Vector{Int}
    basis_indices::SVector{N,Int}
    L::SMatrix{N,N,Int}
    Q::Symmetric{Float64}
end

"""
        SaddlePoint

Represents the saddle point formed at the moment two drainage basins meet.

Fields:
- `labels`: tuple of basin labels assigned during the watershed pre-processing step.
    The labels are stored in ascending order: `labels[1] <= labels[2]`.
- `sites`: indices of the neighboring sites where the basins meet.
- `values`: density values at these sites.
"""
struct SaddlePoint
    labels::Tuple{Int,Int}
    sites::Tuple{Int,Int}
    values::Tuple{Float64,Float64}

    function SaddlePoint(
        labels::Tuple{Int,Int},
        sites::Tuple{Int,Int},
        values::Tuple{Float64,Float64},
    )
        if labels[1] <= labels[2]
            return new(labels, sites, values)
        end
        return new((labels[2], labels[1]), (sites[2], sites[1]), (values[2], values[1]))
    end
end


"""
        WatershedResult{N}

Represents the result of the watershed segmentation algorithm.

Fields:
- `wg`: the underlying watershed grid parameters.
- `ρ`: the values of the density at the grid sites.
- `labels`: the labels assigned to grid sites on the pre-processing pass.
- `summits`: the indices of the highest point for each label.
- `saddles`: a dictionary of the saddle points indexed by ordered tuple of labels.
- `basins`: the indices of drainage basins, constructed on the postprocessing stage.
    Each basin may comprise sites with several labels; the basin index equals
    the smallest label of the sites composing the basin.
"""
struct WatershedResult{N}
    wg::WatershedGrid{N}
    ρ::Vector{Float64}
    labels::Vector{Int}
    summits::Vector{Int}
    saddles::Dict{Tuple{Int,Int}, SaddlePoint}
    basins::Vector{Int}
end



"""
    WatershedResult(wg::WatershedGrid{N}) where N

Construct an empty `WatershedResult{N}` from a `WatershedGrid`.

Allocates storage for the density `values` and `labels` arrays (sized to the
grid), and initializes `summits`, `saddles`, and `basins` as empty collections.
The fields are intended to be filled in by the watershed algorithm.
"""
function WatershedResult(wg::WatershedGrid{N}) where N
    n_sites = wg.grid.size[1]
    return WatershedResult{N}(
        wg,
        Vector{Float64}(undef, n_sites),
        zeros(Int, n_sites)
        Int[],
        Dict{Tuple{Int,Int}, SaddlePoint}(),
        Int[]
    )
end

