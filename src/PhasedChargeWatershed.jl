module PhasedChargeWatershed

import StaticArrays: SMatrix, SVector, MMatrix
import LinearAlgebra: norm, transpose, inv, ⋅, Symmetric, cholesky, svd, qr, mul!, diag
import FFTW: irfft, plan_irfft
import NormalForms: snf
import ChargeFlipPhaser: basis_of_dense_packing
import SpaceGroups: SpaceGroupQuotient, make_orbit, ExtinctOrbit, ComplexOrbit, RealOrbit, @SGE
import JSON: parsefile
import JLD2: jldsave, load
import Makie: Figure, Axis, DataAspect, reset_limits!, hidedecorations!, GridLayout, Label,
    Textbox, Slider, Button, Menu, set_close_to!, Observable, lift, on, notify, heatmap!,
    Colorbar, colsize!, Aspect, resize_to_layout!

export PhasedPeak, PhasedData, SamplingGrid, WatershedGrid, SaddlePoint, WatershedResult
export load_data, save_result, load_result
export sample_density, pre_watershed, sample_pre_watershed_labels, update_basins!
export build_viewer, build_basin_controls, global_density_limits

include("types.jl")
include("nearest_neighbors.jl")
include("io.jl")
include("sampling.jl")
include("watershed.jl")
include("gui.jl")

end