@testset "Kaiser normalization" begin
    m = Matrix{Float64}(I(3))
    @test kaiser_normalize(m)[1] == I(3)
    @test kaiser_normalize(m)[2] == ones(3)
    @test kaiser_denormalize(m, fill(0.5, 3)) == 0.5 * m

    m = copy(A)
    _, w = kaiser_normalize!(m)
    @test all(w .<= 1)
    @test all(norm.(eachrow(m)) .≈ 1)

    @test kaiser_denormalize!(m, w) ≈ A
end
