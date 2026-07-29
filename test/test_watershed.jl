using PhasedChargeWatershed: load_data, pre_watershed, sample_pre_watershed_labels,
    update_basins!

@testset "Watershed tests" begin
    pd = load_data(joinpath(@__DIR__, "data", "synthetic2.json"))
    result = pre_watershed(pd; n_attempts = 100)
    average(a, b) = (a + b) / 2

    # sampled at the watershed sites themselves, the labels must come back unchanged
    @test sample_pre_watershed_labels(result, result.wg.grid) == result.labels

    # with both controls at rest the pre-watershed labelling is reproduced
    update_basins!(result, 0, 1, average)
    @test result.basins == collect(eachindex(result.summits))

    # only the highest summit reaches a threshold of 100%
    update_basins!(result, 100, 1, average)
    @test count(!iszero, result.basins) == 1

    # a fused group is indexed by its smallest member
    update_basins!(result, 0, 0, average)
    groups = Dict{Int,Vector{Int}}()
    for (i, b) in enumerate(result.basins)
        b == 0 || push!(get!(groups, b, Int[]), i)
    end
    @test all(k == minimum(v) for (k, v) in groups)
end
