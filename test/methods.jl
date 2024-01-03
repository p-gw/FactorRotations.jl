function test_criterion_and_gradient(method, Λ)
    Q, ∇Q = criterion_and_gradient(method, Λ)

    @test Q isa Real
    @test ∇Q isa AbstractMatrix{<:Real}
    @test size(∇Q) == size(Λ)
    return nothing
end

@testset "factor rotation methods" begin
    @testset "Biquartimax" begin
        method = Biquartimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Biquartimax is a special case of Oblimin
        @test rotate(A, Biquartimax()) ≈ rotate(A, Oblimin(gamma = 0.5, orthogonal = true))
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

        @test rotate(A, CrawfordFerguson(kappa = 0, orthogonal = true)) ≈
              rotate(A, Quartimax())

        @test rotate(A, CrawfordFerguson(kappa = 0, orthogonal = true)) ≈
              rotate(A, Oblimin(gamma = 0, orthogonal = true))

        @test rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true)) ≈
              rotate(A, Varimax())

        @test rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true)) ≈
              rotate(A, Oblimin(gamma = 1, orthogonal = true))

        # TODO: Equamax: kappa = k/2p
        # TODO: Parsimax: kappa = (k - 1)/(p + k - 2)
        # TODO: Factor Parsimony: kappa = 1

        # oblique case
        method = CrawfordFerguson(kappa = 0.5, orthogonal = false)
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

    @testset "Oblimin" begin
        # orthogonal case
        method = Oblimin(gamma = 0.5, orthogonal = true)
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # oblique case
        method = Oblimin(gamma = 0.0, orthogonal = false)
        @test isoblique(method)
    end

    @testset "Quartimax" begin
        method = Quartimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Quartimax is a special case of Oblimin
        @test rotate(A, Quartimax()) ≈ rotate(A, Oblimin(gamma = 0, orthogonal = true))
    end

    @testset "Varimax" begin
        method = Varimax()
        @test isorthogonal(method)
        test_criterion_and_gradient(method, A)

        # Varimax is a special case of Oblimin
        @test rotate(A, Varimax()) ≈ rotate(A, Oblimin(gamma = 1, orthogonal = true))

        # Varimax is a special case of Crawford-Ferguson
        p = size(A, 1)
        @test rotate(A, Varimax()) ≈
              rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true))
    end
end
