function test_criterion_and_gradient(method, Λ)
    Q, ∇Q = criterion_and_gradient(method, Λ)

    @test Q isa Real
    @test ∇Q isa AbstractMatrix{<:Real}
    @test size(∇Q) == size(Λ)
    return nothing
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
        @test rotate(A, Biquartimax(); init) ≈
              rotate(A, Oblimin(gamma = 0.5, orthogonal = true); init)
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
        componentloss = rotate(A, ComponentLoss(x -> x^4, orthogonal = true))
        quartimax = rotate(A, Quartimax())
        @test isapprox(componentloss, quartimax, atol = 1e-5)

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

        @test rotate(A, CrawfordFerguson(kappa = 0, orthogonal = true); init) ≈
              rotate(A, Quartimax(); init)

        @test rotate(A, CrawfordFerguson(kappa = 0, orthogonal = true); init) ≈
              rotate(A, Oblimin(gamma = 0, orthogonal = true); init)

        @test rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true); init) ≈
              rotate(A, Varimax(); init)

        @test rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true); init) ≈
              rotate(A, Oblimin(gamma = 1, orthogonal = true); init)

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
        oblimax = rotate(A, method; init)
        quartimax = rotate(A, Quartimax(); init)
        @test isapprox(oblimax, quartimax, atol = 1e-6)

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
        Ar = rotate(A, Quartimax(); init)
        @test Ar ≈ rotate(A, Oblimin(gamma = 0, orthogonal = true); init)

        # test that rotation result is identical published results by
        # Bernaards & Jennrich (2005) within the reported accuracy of 7 digits
        Ar = rotate(A, Quartimax(); init, atol = 1e-5)
        Ar = round.(Ar, digits = 7)

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
        @test rotate(A, Varimax(); init) ≈
              rotate(A, Oblimin(gamma = 1, orthogonal = true); init)

        # Varimax is a special case of Crawford-Ferguson
        p = size(A, 1)
        @test rotate(A, Varimax(); init) ≈
              rotate(A, CrawfordFerguson(kappa = 1 / p, orthogonal = true); init)
    end
end
