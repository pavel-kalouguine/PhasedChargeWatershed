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
    gram = round.(Int, B' * B)
    R = cholesky(Symmetric(Float64.(gram))).U
    norm_bound = minimum(gram[i, i] for i in 1:n)
    tol = 1e-6

    found = SVector{n,Int}[]
    x = zeros(Int, n)
    search!(found, x, n, 0.0, R, gram, norm_bound, tol)
    isempty(found) && return found
    q_min = minimum(v' * gram * v for v in found)
    filter(v -> v' * gram * v == q_min, found)
end
