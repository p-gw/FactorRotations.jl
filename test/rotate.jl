@testset "rotate" begin
    # initial values
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(10, 10))
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(8, 2))

    @test rotate(A, Varimax()) == rotate(A, Varimax(); init)

    # in-place rotation
    B = copy(A)
    rotate!(B, Quartimax())
    @test B == rotate(A, Quartimax())
end
