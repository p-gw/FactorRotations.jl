@testset "utility functions" begin
    @testset "zerodiag" begin
        m = ones(3, 3)
        @test FactorRotations.zerodiag!(m) == ones(3, 3) - I
        @test diag(m) == zeros(3)
    end

    @testset "nthsmallest" begin
        @test FactorRotations.nthsmallest(1:10, 6) == 6
        @test FactorRotations.nthsmallest(1:10, 1) == 1
        @test FactorRotations.nthsmallest([1, 10, 100], 3) == 100
        @test FactorRotations.nthsmallest([1 3; 7 13], 2) == 3
    end

    @testset "random_orthogonal_matrix" begin
        m = FactorRotations.random_orthogonal_matrix(10)
        @test size(m) == (10, 10)
    end

    @testset "setverbosity!" begin
        @test FactorRotations.VERBOSITY[] == false
        @test_logs (:info, "Logging is disabled globally.") setverbosity!(false)
        @test_logs (:info, "Logging is enabled globally.") setverbosity!(true)
        @test FactorRotations.VERBOSITY[] == true

        setverbosity!(false)  # disable logging for following tests
    end
end
