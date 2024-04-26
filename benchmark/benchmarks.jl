using BenchmarkTools
using FactorRotations

const SUITE = BenchmarkGroup()

SUITE["criterion_and_gradient!"] = BenchmarkGroup()

methods = [
    Biquartimax(),
    Biquartimin(),
    CrawfordFerguson(kappa = 0.5),
    Geomin(),
    Infomax(),
    MinimumEntropyRatio(),
    MinimumEntropy(),
    Oblimax(),
    Oblimin(gamma = 0.5),
    Quartimax(),
    Simplimax(m = 5),
    TandemCriterionI(),
    TandemCriterionII(),
    # TargetRotation(zeros(10, 3)),  exclude due to printing issues
    Varimax(),
]

for method in methods
    SUITE["criterion_and_gradient!"][method] = @benchmarkable(
        criterion_and_gradient!($(zeros(10, 3)), $method, $(rand(10, 3))),
        evals = 10,
        samples = 1000,
    )
end
