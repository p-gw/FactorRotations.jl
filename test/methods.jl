function test_criterion_and_gradient(method, Λ)
    Q, ∇Q = criterion_and_gradient(method, Λ)

    @test Q isa Real
    @test ∇Q isa AbstractMatrix{<:Real}
    @test size(∇Q) == size(Λ)
    return nothing
end

function test_rotate(method, Λ; init)
    rot = rotate(Λ, method; init)
    p, k = size(Λ)

    @test size(loadings(rot)) == (p, k)
    @test size(rotation(rot)) == (k, k)
    @test size(factor_correlation(rot)) == (k, k)
    @test loadings(rot) * rotation(rot)' ≈ Λ

    if isorthogonal(method)
        @test factor_correlation(rot) ≈ I
    end
end

function test_equivalence(Λ, m1::RotationMethod, m2::RotationMethod; kwargs...)
    r1 = rotate(Λ, m1; kwargs...)
    r2 = rotate(Λ, m2; kwargs...)

    @test isapprox(loadings(r1), loadings(r2), atol = 1e-5)
    @test isapprox(rotation(r1), rotation(r2), atol = 1e-5)
    @test isapprox(factor_correlation(r1), factor_correlation(r2), atol = 1e-5)
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

        method = Absolmin(0)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
        @test_throws ArgumentError Absolmin(-1.0)
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

        # TODO: Equamax: kappa = k/2p
        # TODO: Parsimax: kappa = (k - 1)/(p + k - 2)
        # TODO: Factor Parsimony: kappa = 1

        # oblique case
        method = CrawfordFerguson(kappa = 0.5, orthogonal = false)
        @test isoblique(method)
        test_criterion_and_gradient(method, A)
    end

    @testset "Cubimax" begin
        method = Cubimax()
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
        # orthogonal case
        method = Oblimin(gamma = 0.5, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique case
        method = Oblimin(gamma = 0.0, orthogonal = false)
        @test isoblique(method)
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
        @test_throws ArgumentError criterion_and_gradient(TargetRotation([0 1; 1 0]), A)

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

        # test that rotation result is identical published results by
        # Bernaards & Jennrich (2005) within the reported accuracy of 7 digits
        Ar = rotate(A, Quartimax(); init, atol = 1e-5)
        Ar = round.(loadings(Ar), digits = 7)

        pub = [
            0.8987554 0.1948197
            0.9339440 0.1297446
            0.9021319 0.1038604
            0.8765090 0.1712805
            0.3155758 0.8764747
            0.2511265 0.7734879
            0.1980102 0.7146775
            0.3078601 0.6593331
        ]

        @test Ar ≈ pub
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
end
