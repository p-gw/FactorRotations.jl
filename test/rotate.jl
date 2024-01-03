@testset "rotate" begin
    # initial values
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(10, 10))
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(8, 2))
end
