# Command-line script running the pre-watershed segmentation.
#
#   julia --project prewatershed.jl <input.json> <output.jld2> [density_factor] [n_attempts]
#
# Loads the phased data from <input.json>, runs the pre-watershed segmentation on it and
# saves the resulting `WatershedResult` into <output.jld2> under the key "result". The
# result carries the phased data along with it, so the saved file is self-contained: the
# density and the basins can never end up coming from two unrelated files.
#
# Both numeric arguments are optional and are simply forwarded to `pre_watershed`, which
# validates them and uses the same defaults: `density_factor` (1.0) controls how dense the
# watershed grid is, `n_attempts` (1000) how many candidate lattices are tried.
# The segmentation is the slow part, so the script is meant to be run once per set of
# parameters, the result being reused later.
#
# The saved result is read back with
#   using JLD2
#   result = load("output.jld2", "result")

using JLD2
import PhasedChargeWatershed: load_data, pre_watershed

if !(2 <= length(ARGS) <= 4)
    error("usage: julia --project prewatershed.jl <input.json> <output.jld2> [density_factor] [n_attempts]")
end
input_path, output_path = ARGS[1], ARGS[2]
isfile(input_path) || error("input file not found: $input_path")
output_dir = dirname(abspath(output_path))
isdir(output_dir) || error("output directory does not exist: $output_dir")

# Read a trailing optional argument, falling back to the default
function read_optional(position, T, name, default)
    length(ARGS) < position && return default
    value = tryparse(T, ARGS[position])
    value === nothing && error("$name: expected $T, got \"$(ARGS[position])\"")
    return value
end

density_factor = read_optional(3, Float64, "density_factor", 1.0)
n_attempts = read_optional(4, Int, "n_attempts", 1000)

phased_data = load_data(input_path)
println("loaded $(length(phased_data.peaks)) peaks from $input_path")

println("running the pre-watershed with density_factor = $density_factor, n_attempts = $n_attempts ...")
result = pre_watershed(phased_data; density_factor = density_factor, n_attempts = n_attempts)
println("got $(length(result.summits)) basins and $(length(result.saddles)) saddle points",
        " on a grid of $(result.wg.grid.size[1]) sites")

jldsave(output_path; result = result)
println("saved the result to $output_path")
