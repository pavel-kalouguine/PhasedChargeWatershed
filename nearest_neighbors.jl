using StaticArrays, LinearAlgebra
import ChargeFlipPhaser: basis_of_dense_packing

function nearest_neighbors(n::Int)
    B = basis_of_dense_packing(n)
    G = B' * B
    gram = round.(Int, G)
    gram_inv = inv(G)
    min_diag = minimum(gram[i, i] for i in 1:n)
    bounds = [floor(Int, sqrt(min_diag * gram_inv[i, i])) + 1 for i in 1:n]
    ranges = [-b:b for b in bounds]
    best = typemax(Int)
    neighbors = SVector{n,Int}[]
    for t in Iterators.product(ranges...)
        x = SVector{n,Int}(t)
        iszero(x) && continue
        q = x' * gram * x
        if q < best
            best = q
            empty!(neighbors)
            push!(neighbors, x)
        elseif q == best
            push!(neighbors, x)
        end
    end
    neighbors
end
