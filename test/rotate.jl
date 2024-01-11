@testset "rotate" begin
    # initial values
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(10, 10))
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(8, 2))

    @test rotate(A, Varimax()) == rotate(A, Varimax(); init)

    # in-place rotation
    B = copy(A)
    rotate!(B, Quartimax())
    @test B == rotate(A, Quartimax())

    @testset "RotationState" begin
        orthogonal_state = FactorRotations.RotationState(Orthogonal, init, A)
        @test orthogonal_state.init == init
        @test orthogonal_state.A == A
        @test orthogonal_state.T == init
        @test orthogonal_state.L == A * init
        @test orthogonal_state.iterations == FactorRotations.IterationState[]

        oblique_state = FactorRotations.RotationState(Oblique, init, A)
        @test oblique_state.init == init
        @test oblique_state.A == A
        @test oblique_state.T == init
        @test oblique_state.Ti == inv(init)
        @test oblique_state.L == A * inv(init)'
        @test oblique_state.iterations == FactorRotations.IterationState[]
    end

    @testset "gradient_f" begin
        orthogonal_state = FactorRotations.RotationState(Orthogonal, init, A)
        @test FactorRotations.gradient_f(orthogonal_state, zeros(size(A))) ==
              zeros(size(init))

        oblique_state = FactorRotations.RotationState(Oblique, init, A)
        @test FactorRotations.gradient_f(oblique_state, zeros(size(A))) == zeros(size(init))
    end

    @testset "project_G!" begin
        Gp = zeros(Float64, size(init))
        G = rand(size(init)...)
        orthogonal_state = FactorRotations.RotationState(Orthogonal, init, A)
        @test FactorRotations.project_G!(orthogonal_state, Gp, G) != zeros(size(init))
        @test Gp != zeros(size(init))

        Gp = zeros(Float64, size(init))
        oblique_state = FactorRotations.RotationState(Oblique, init, A)
        @test FactorRotations.project_G!(oblique_state, Gp, G) != zeros(size(init))
        @test Gp != zeros(size(init))
    end

    @testset "update_state!" begin
        Tt = [1 1; 1 0]

        orthogonal_state = FactorRotations.RotationState(Orthogonal, init, A)
        @test FactorRotations.update_state!(orthogonal_state, Tt) == A * Tt
        @test orthogonal_state.L == A * Tt

        oblique_state = FactorRotations.RotationState(Oblique, init, A)
        @test FactorRotations.update_state!(oblique_state, Tt) == A * inv(Tt)'
        @test oblique_state.Ti == inv(Tt)
        @test oblique_state.L == A * inv(Tt)'
    end
end
