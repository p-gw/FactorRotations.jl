using FactorRotations
using LinearAlgebra
using Test

A = [
    0.830 -0.396
    0.818 -0.469
    0.777 -0.470
    0.798 -0.401
    0.786 0.500
    0.672 0.458
    0.594 0.444
    0.647 0.333
];

init = Matrix{Float64}(I, 2, 2)

@testset "FactorRotations.jl" begin
    include("methods.jl")
    include("rotate.jl")
end
