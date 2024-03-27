@testset "reflect" begin
    @test FactorRotations.reflect_cols(ones(8, 3)) == ones(1, 3)
    @test FactorRotations.reflect_cols(hcat(ones(8), -ones(8))) == [1 -1]
    @test FactorRotations.reflect_cols(hcat(-ones(8), zeros(8), ones(8))) == [-1 1 1]

    @testset "FactorRotation" begin
        Imat = diagm(ones(Float64, 3))
        r = FactorRotation(hcat(ones(4), zeros(4), -ones(4)), Imat, ones(4))
        expectation = FactorRotation(
            hcat(ones(4), zeros(4), ones(4)),
            diagm(Float64[1, 1, -1]),
            ones(4),
        )

        @test reflect(r).T == expectation.T
        @test reflect(r).L == expectation.L
        @test reflect(r).phi == expectation.phi

        reflect!(r)
        @test r.T == expectation.T
        @test r.L == expectation.L
        @test r.phi == expectation.phi
    end
end
