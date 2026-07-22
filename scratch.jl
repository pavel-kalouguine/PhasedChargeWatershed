using StaticArrays

include(joinpath(@__DIR__, "src", "PhasedChargeWatershed.jl"))

import PhasedChargeWatershed: load_data, SamplingGrid, sample_density, PhasedData, PhasedPeak, nearest_neighbors
import ChargeFlipPhaser: basis_of_dense_packing
import NormalForms: snf
import LinearAlgebra: svd, inv, transpose, norm, Symmetric, det

function quadratic_moment(peaks::Vector{PhasedPeak{N}}) where N
    # Compute the quadratic moment of the phased peaks
    Q = zeros(Float64, N, N)
    for peak in peaks
        k = peak.k
        f = peak.f
        Q += abs2(f) * (k * k')
    end
    return Symmetric(Q) / sum(abs2.(getfield.(peaks, :f)))
end



filepath = joinpath(@__DIR__, "data", "synthetic.json")
phased_data = load_data(filepath)


function create_watershed_grid(phased_data::PhasedData{N,D}) where {N,D}
    Q = quadratic_moment(phased_data.peaks)
    B0=basis_of_dense_packing(N)
    neighbors = nearest_neighbors(B0)
    scaling_factor=10.0 # Linear scaling
    n_attempts=1
    i=0
    while true
        i+=1
        R = svd(randn(N, N)).U # Random orthogonal matrix
        B=R*B0 # Rotated basis
        L=round.(Int, scaling_factor*inv(B)*Q^0.5)
        #L=scaling_factor*inv(B)*Q^0.5
        for u in neighbors
            Li=round.(Int, inv(L)*det(L)).//round(Int, det(L))
            #v=inv(L)*u
            v=Li*u
            #println(v'*Q*v)
            println(v)
        end
        println(L)
        if i>=n_attempts
            break
        end
        
    end
end

create_watershed_grid(phased_data)