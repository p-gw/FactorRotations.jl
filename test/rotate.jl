@testset "rotate" begin
    # initial values
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(10, 10))
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(8, 2))

    rot_default_init = rotate(A, Varimax())
    rot_identity_init = rotate(A, Varimax(); init)
    @test loadings(rot_default_init) ≈ loadings(rot_identity_init)
    @test rotation(rot_default_init) ≈ rotation(rot_identity_init)
    @test factor_correlation(rot_default_init) ≈ factor_correlation(rot_identity_init)

    # in-place rotation
    B = copy(A)
    rotate!(B, Quartimax())
    @test B == loadings(rotate(A, Quartimax()))

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

    @testset "project_X" begin
        state = FactorRotations.RotationState(Orthogonal, init, A)
        @test FactorRotations.project_X(state, I(2)) == I(2)

        state = FactorRotations.RotationState(Oblique, init, A)
        @test FactorRotations.project_X(state, I(2)) == I(2)
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

    @testset "MultivariateStatsExt" begin
        using MultivariateStats

        X = rand(3, 100)
        methods = [FactorAnalysis, PCA, PPCA]

        for method in methods
            model = fit(FactorAnalysis, X)
            raw_loadings = MultivariateStats.loadings(model)
            rot = rotate(raw_loadings, Varimax())
            rotated_loadings = FactorRotations.loadings(rot)

            @test size(raw_loadings) == (3, 2)
            @test size(FactorRotations.loadings(rotate(model, Varimax()))) == (3, 2)
            @test rotated_loadings == FactorRotations.loadings(rotate(model, Varimax()))
            @test MultivariateStats.loadings(model) == raw_loadings

            @test rotated_loadings == rotate!(model, Varimax())
            @test MultivariateStats.loadings(model) == rotated_loadings
        end
    end
end
