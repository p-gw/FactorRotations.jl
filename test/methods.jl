function test_criterion_and_gradient(method, Λ)
    ∇Q = fill!(similar(Λ), NaN)
    Q = @inferred(criterion_and_gradient!(∇Q, method, Λ))
    @test Q isa Real
    @test all(isfinite, ∇Q)

    # test criterion-only calculation
    Q2 = @inferred(criterion_and_gradient!(nothing, method, Λ))
    @test Q2 == Q

    # test criterion() wrapper
    Q3 = @inferred(criterion(method, Λ))
    @test Q3 isa Real
    @test Q3 == Q

    return nothing
end

function test_rotate(method, Λ; init)
    rot = @inferred(rotate(Λ, method; init))
    p, k = size(Λ)

    @test size(@inferred(loadings(rot))) == (p, k)
    @test size(@inferred(rotation(rot))) == (k, k)
    @test size(@inferred(factor_correlation(rot))) == (k, k)
    @test loadings(rot) * rotation(rot)' ≈ Λ

    if isorthogonal(method)
        @test factor_correlation(rot) ≈ I
    end
end

function test_equivalence(Λ, m1::RotationMethod, m2::RotationMethod; kwargs...)
    r1 = rotate(Λ, m1; kwargs...)
    r2 = rotate(Λ, m2; kwargs...)

    @test loadings(r1) ≈ loadings(r2) atol = 1e-5
    @test rotation(r1) ≈ rotation(r2) atol = 1e-5
    @test factor_correlation(r1) ≈ factor_correlation(r2) atol = 1e-5
end

@testset "factor rotation autodiff fallback" begin
    method = ComponentLoss(abs2, orthogonal = true)
    ∇Q = fill!(similar(A), NaN)

    FactorRotations.set_autodiff_backend(:ABC)
    @test_throws "ABC autodiff backend is not supported" criterion_and_gradient!(∇Q, method, A)
    FactorRotations.set_autodiff_backend(:Enzyme)
    #@test_throws "Enzyme.jl autodiff backend is not loaded" criterion_and_gradient!(∇Q, method, A)
end

@testset "factor rotation methods" begin
    @testset "utility functions" begin
        @test isorthogonal(Varimax()) != isoblique(Varimax())
        @test isorthogonal(Oblimax(orthogonal = false)) !=
              isoblique(Oblimax(orthogonal = false))
    end

    @testset "Biquartimax" begin
        method = Biquartimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Biquartimax is a special case of Oblimin
        test_equivalence(A, Biquartimax(), Oblimin(gamma = 0.5, orthogonal = true); init)
    end

    @testset "Biquartimin" begin
        # Example 3.1 in Jennrich & Bentler (2011)
        bifactorA = [
            1.17 0.78 0.18
            2.08 0.78 -0.22
            1.17 0.78 0.18
            2.15 -0.62 -0.08
            1.23 -0.62 0.32
            2.15 -0.62 -0.08
        ]

        # oblique case
        method = Biquartimin()
        @test isoblique(method)
        test_criterion_and_gradient(method, bifactorA)

        # orthogonal case
        method = Biquartimin(orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, bifactorA)
    end

    @testset "ComponentLoss" begin
        # orthogonal case
        method = ComponentLoss(abs2, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        method = KatzRohlf(0.3)
        @test isorthogonal(method)
        @test_throws ArgumentError KatzRohlf(-1.0)
        @test_throws ArgumentError KatzRohlf(0)
        test_criterion_and_gradient(method, A)

        method = LinearRightConstant(0.3)
        @test isorthogonal(method)
        @test_throws ArgumentError LinearRightConstant(-1.0)
        @test_throws ArgumentError LinearRightConstant(0)
        test_criterion_and_gradient(method, A)

        # component loss identical to quartimax
        test_equivalence(A, ComponentLoss(x -> x^4, orthogonal = true), Quartimax(); init)

        # oblique case
        method = ComponentLoss(abs2, orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)

        method = Concave(1)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
        @test_throws ArgumentError Concave(0)
        @test_throws ArgumentError Concave(-2)

        method = Absolmin(1e-5)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
        @test_throws ArgumentError Absolmin(-1.0)
        @test_throws ArgumentError Absolmin(0)
    end

    @testset "CrawfordFerguson" begin
        @test_throws ArgumentError CrawfordFerguson(kappa = 2.0)
        @test_throws ArgumentError CrawfordFerguson(kappa = -0.2)

        # orthogonal case
        method = CrawfordFerguson(kappa = 0.2, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # for the orthogonal case, Crawford-Ferguson and Oblimin are equivalent
        p, k = size(A)

        test_equivalence(
            A,
            CrawfordFerguson(kappa = 0, orthogonal = true),
            Quartimax();
            init,
        )
        test_equivalence(
            A,
            CrawfordFerguson(kappa = 0, orthogonal = true),
            Oblimin(gamma = 0, orthogonal = true);
            init,
        )
        test_equivalence(
            A,
            CrawfordFerguson(kappa = 1 / p, orthogonal = true),
            Varimax();
            init,
        )
        test_equivalence(
            A,
            CrawfordFerguson(kappa = 1 / p, orthogonal = true),
            Oblimin(gamma = 1, orthogonal = true);
            init,
        )
        test_equivalence(
            A,
            Equamax(),
            CrawfordFerguson(kappa = k / (2p), orthogonal = true);
            init,
        )
        test_equivalence(
            A,
            Parsimax(),
            CrawfordFerguson(kappa = (k - 1) / (p + k - 2), orthogonal = true);
            init,
        )

        # TODO: Factor Parsimony: kappa = 1

        # oblique case
        method = CrawfordFerguson(kappa = 0.5, orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Equamax" begin
        method = Equamax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Geomin" begin
        @test_throws ArgumentError Geomin(epsilon = -1)

        method = Geomin()
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Infomax" begin
        # orthogonal case
        method = Infomax(orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique case
        method = Infomax(orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "MinimumEntropy" begin
        method = MinimumEntropy()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "MinimumEntropyRatio" begin
        method = MinimumEntropyRatio()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Oblimax" begin
        # orthogonal case
        method = Oblimax(orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Oblimax is equivalent to Quartimax in the orthogonal case
        test_equivalence(A, method, Quartimax(); init)

        # oblique case
        method = Oblimax(orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Oblimin" begin
        p, k = size(A)

        # orthogonal case
        method = Oblimin(gamma = 0.5, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique case
        method = Oblimin(gamma = 0.0, orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)

        test_equivalence(A, Equamax(), Oblimin(gamma = k / 2, orthogonal = true); init)
        test_equivalence(
            A,
            Parsimax(),
            Oblimin(gamma = p * (k - 1) / (p + k - 2), orthogonal = true);
            init,
        )
    end

    @testset "Parsimax" begin
        method = Parsimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "PatternSimplicity" begin
        # orthogonal case
        method = PatternSimplicity(orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique case
        method = PatternSimplicity(orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "TandemCriteria" begin
        method = TandemCriterionI()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        method = TandemCriterionII()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        @test_throws ArgumentError TandemCriteria(keep = 0)
        @test_throws ArgumentError TandemCriteria(keep = -1)
        @test_throws ArgumentError TandemCriteria(keep = 1)

        method = TandemCriteria(keep = 2)
        @test isorthogonal(method)
        rot = rotate(A, method)
        @test size(loadings(rot)) == (8, 2)
    end

    @testset "TargetRotation" begin
        @test_throws ArgumentError criterion_and_gradient!(
            similar(A),
            TargetRotation([0 1; 1 0]),
            A,
        )

        # orthogonal + complete case
        method = TargetRotation(similar(A), orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # orthogonal + partial specification
        target = [missing 1; 0 1; 0 1; 0 1; 1 0; 1 0; 1 0; 1 0]
        method = TargetRotation(target, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique + complete case
        method = TargetRotation(similar(A), orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)

        # oblique + partial specification
        method = TargetRotation(target, orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Quartimax" begin
        method = Quartimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Quartimax is a special case of Oblimin
        test_equivalence(A, Quartimax(), Oblimin(gamma = 0, orthogonal = true); init)

        # test that rotation result is identical to GPArotation
        Ar = rotate(A, Quartimax(); init, g_atol = 1e-7)

        # loadings published in Bernaards & Jennrich (2005)
        # within the reported accuracy of 7 digits
        #pub = [
        #    0.8987554 0.1948197
        #    0.9339440 0.1297446
        #    0.9021319 0.1038604
        #    0.8765090 0.1712805
        #    0.3155758 0.8764747
        #    0.2511265 0.7734879
        #    0.1980102 0.7146775
        #    0.3078601 0.6593331
        #]
        # loadings obtained with GPArotation v2024.3
        # quartimax(A, eps=1e-8, maxit=50000, randomStarts=10)
        gpa = [
            0.8987545678868889 0.19482357840480186
            0.9339434064715286 0.1297486551312077
            0.9021314838553653 0.10386426641014213
            0.8765082522883497 0.1712842189765975
            0.31557202157519415 0.87647606881132
            0.2511231928032839 0.7734889411208703
            0.19800711751346906 0.7146783762042948
            0.3078572424280441 0.6593344510069232
        ]

        @test FactorRotations.loadings(Ar) ≈ gpa atol = 1e-6
        @test criterion(Quartimax(), Ar.L) ≈ criterion(Quartimax(), gpa) atol = 1e-8
    end

    @testset "Simplimax" begin
        @test_throws ArgumentError Simplimax(m = 0)
        @test_throws ArgumentError Simplimax(m = -10)

        method = Simplimax(m = 8)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Varimax" begin
        method = Varimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Varimax is a special case of Oblimin
        test_equivalence(A, Varimax(), Oblimin(gamma = 1, orthogonal = true); init)

        # Varimax is a special case of Crawford-Ferguson
        p = size(A, 1)
        test_equivalence(
            A,
            Varimax(),
            CrawfordFerguson(kappa = 1 / p, orthogonal = true);
            init,
        )
    end

    @testset "Missing criterion implementation" begin
        struct NoCriterion <: RotationMethod{Orthogonal} end

        @test_throws "NoCriterion does not implement" criterion(NoCriterion(), randn(6, 6))
        # Enzyme.jl would refuse to autodiff because it detects that fallback criterion() throws an error
        @test_throws "NoCriterion does not implement" criterion_and_gradient!(
            randn(6, 5),
            NoCriterion(),
            randn(6, 5),
        )
    end
end
