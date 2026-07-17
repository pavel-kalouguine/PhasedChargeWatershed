module PhasedChargeWatershed

import StaticArrays: SMatrix, SVector, MMatrix
import LinearAlgebra: norm, transpose, inv, ⋅, Symmetric, cholesky, svd, qr
import FFTW: irfft
import ChargeFlipPhaser: basis_of_dense_packing
import SpaceGroups: SpaceGroupQuotient, make_orbit, ExtinctOrbit, ComplexOrbit, RealOrbit, @SGE
import JSON: parsefile  


include("types.jl")
include("nearest_neighbors.jl")
include("io.jl")
include("sampling.jl")

end