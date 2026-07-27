using Test
using PhasedChargeWatershed

@testset "PhasedChargeWatershed.jl" begin
    include("test_types.jl")
    include("test_io.jl")
    include("test_watershed.jl")
end
