using LinearAlgebra: norm, I

@testset "nearest_neighbors tests" begin

    @testset "cubic lattice (dimensions 1-10)" begin
        for n in 1:10
            B = Matrix{Float64}(2I, n, n)
            neighbors = PhasedChargeWatershed.nearest_neighbors(B)

            # The cubic lattice with spacing 2 has 2n nearest neighbors:
            # one at +2 and one at -2 along each coordinate axis.
            @test length(neighbors) == 2n

            # All of them must be at distance exactly 2 from the origin.
            for v in neighbors
                @test norm(B * v) ≈ 2.0
            end
        end
    end

    @testset "densest known sphere packings (dimensions 1-8)" begin
        # Known kissing numbers of the densest sphere packings in dimensions 1 to 8.
        kissing_numbers = Dict(1 => 2, 2 => 6, 3 => 12, 4 => 24, 5 => 40, 6 => 72, 7 => 126, 8 => 240)

        for n in 1:8
            B = PhasedChargeWatershed.basis_of_dense_packing(n)
            neighbors = PhasedChargeWatershed.nearest_neighbors(B)

            @test length(neighbors) == kissing_numbers[n]

            for v in neighbors
                @test norm(B * v) ≈ 2.0
            end
        end
    end

end
