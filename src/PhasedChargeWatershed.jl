module PhasedChargeWatershed

import StaticArrays: SMatrix, SVector, MMatrix
import LinearAlgebra: norm, transpose, inv, ⋅, Symmetric, cholesky, svd, qr
import FFTW: irfft
import ChargeFlipPhaser: basis_of_dense_packing
import SpaceGroups: SpaceGroupQuotient, make_orbit, ExtinctOrbit, ComplexOrbit, RealOrbit, @SGE
import JSON: parsefile
import Makie: Figure, Axis, AxisAspect, hidedecorations!, GridLayout, Label,
    Textbox, Slider, Button, set_close_to!, Observable, lift, on, heatmap!


include("types.jl")
include("nearest_neighbors.jl")
include("io.jl")
include("sampling.jl")
include("gui.jl")

end