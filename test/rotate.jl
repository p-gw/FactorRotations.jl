@testset "rotate" begin
    # initial values
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(10, 10))
    @test_throws ArgumentError rotate(A, Varimax(), init = rand(8, 2))

    rot_default_init = rotate(A, Varimax(), g_atol=1e-7)
    rot_identity_init = rotate(A, Varimax(); g_atol=1e-7, init)
    @test loadings(rot_default_init) ≈ loadings(rot_identity_init) atol=1e-7
    @test rotation(rot_default_init) ≈ rotation(rot_identity_init) atol=1e-7
    @test factor_correlation(rot_default_init) ≈ factor_correlation(rot_identity_init) atol=1e-7

    # in-place rotation
    B = copy(A)
    rotate!(B, Quartimax())
    @test B == loadings(rotate(A, Quartimax()))

    # random starts
    @test_warn "Ignoring initial starting values" rotate(
        A,
        Varimax(),
        randomstarts = 3,
        init = I(2),
    )

    @test_nowarn rotate(A, Varimax(), randomstarts = 3)

    # convergence
    struct NonConverging <: RotationMethod{Orthogonal} end

    function FactorRotations.criterion_and_gradient(::NonConverging, m::AbstractMatrix)
        return (1.0, ones(size(m)))
    end

    @test_throws ConvergenceError rotate(A, NonConverging())
    @test_throws ConvergenceError rotate(A, NonConverging(), randomstarts = 3)

    @test_warn "did not converge" rotate(
        ones(8, 2),
        LinearRightConstant(1.0),
        randomstarts = true,
    )

    # normalize
    B = copy(A)
    rotate(B, Varimax(), normalize = true)
    @test B == A

    @testset "parse_randomstarts" begin
        @test FactorRotations.parse_randomstarts(true) == 100
        @test FactorRotations.parse_randomstarts(true; default = 10) == 10
        @test FactorRotations.parse_randomstarts(false) == 0
        @test FactorRotations.parse_randomstarts(9) == 9
        @test_throws ArgumentError FactorRotations.parse_randomstarts(-1)
        @test_throws ArgumentError FactorRotations.parse_randomstarts(0)
    end

    @testset "RotationState" begin
        orthogonal_state = FactorRotations.RotationState(Orthogonal, init, A)
        @test orthogonal_state.init == init
        @test orthogonal_state.A == A
        @test orthogonal_state.T == init
        @test orthogonal_state.L == A * init
        @test orthogonal_state.iterations == FactorRotations.IterationState[]
        @test isnan(FactorRotations.minimumQ(orthogonal_state))

        oblique_state = FactorRotations.RotationState(Oblique, init, A)
        @test oblique_state.init == init
        @test oblique_state.A == A
        @test oblique_state.T == init
        @test oblique_state.Ti == inv(init)
        @test oblique_state.L == A * inv(init)'
        @test oblique_state.iterations == FactorRotations.IterationState[]
        @test isnan(FactorRotations.minimumQ(oblique_state))

        struct BadRotation <: FactorRotations.RotationType end
        @test_throws "Unsupported rotation type BadRotation" FactorRotations.RotationState(BadRotation, init, A)
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
