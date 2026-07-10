using StaticArrays, LinearAlgebra

function search!(found, x, i, acc, R, gram, norm_bound, tol)
    n = length(x)
    if i == 0
        q = x' * gram * x
        (0 < q <= norm_bound) && push!(found, SVector{n,Int}(x))
        return
    end
    offset = sum(R[i, j] * x[j] for j in (i+1):n; init=0.0)
    remaining = norm_bound - acc
    remaining < -tol && return
    radius = sqrt(max(remaining, 0.0) + tol)
    for xi in ceil(Int, (-radius - offset) / R[i, i]):floor(Int, (radius - offset) / R[i, i])
        x[i] = xi
        term = R[i, i] * xi + offset
        search!(found, x, i - 1, acc + term^2, R, gram, norm_bound, tol)
    end
end

function nearest_neighbors(B::AbstractMatrix)
    n = size(B, 1)
    if size(B, 2) != n
        error("B must be a square matrix")
    end
    G = B' * B
    gram = round.(Int, G)
    if maximum(abs.(G - gram)) > 1e-6
        error("B' * B must be an integer matrix")
    end
    R = cholesky(Symmetric(Float64.(gram))).U

    shortest_distance = 2
    norm_bound = shortest_distance^2
    tol = 1e-6
    found = SVector{n,Int}[]
    x = zeros(Int, n)
    search!(found, x, n, 0.0, R, gram, norm_bound, tol)
    if isempty(found)
        error("no lattice vectors found within the specified norm bound")
    end
    q_min = minimum(v' * gram * v for v in found)
    if q_min != norm_bound
        error("found a lattice vector shorter than distance 2")
    end
    filter(v -> v' * gram * v == q_min, found)
end
