using BenchmarkTools
using FactorRotations

const SUITE = BenchmarkGroup()

SUITE["criterion_and_gradient!"] = BenchmarkGroup()

methods = [Varimax(), Quartimax(), Geomin()]

for method in methods
    SUITE["criterion_and_gradient!"][method] = @benchmarkable(
        criterion_and_gradient!($(zeros(10, 3)), $method, $(rand(10, 3))),
        evals = 10,
        samples = 1000,
    )
end
