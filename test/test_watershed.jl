using PhasedChargeWatershed: load_data, pre_watershed, sample_pre_watershed_labels,
    update_basins!

@testset "Watershed tests" begin
    pd = load_data(joinpath(@__DIR__, "data", "synthetic2.json"))
    result = pre_watershed(pd; n_attempts = 100)
    average(a, b) = (a + b) / 2

    # sampled at the watershed sites themselves, the labels must come back unchanged
    @test sample_pre_watershed_labels(result, result.wg.grid) == result.labels

    # a saddle point never rises above the summits it separates, so nothing fuses at 1
    update_basins!(result, 0, 1, average)
    @test all(b == 0 || b == i for (i, b) in enumerate(result.basins))

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
