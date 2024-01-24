# Working with MultivariateStats.jl

FactorRotations.jl provides direct support for models fitted by [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl). Specifically, you can fit a [factor analysis](https://juliastats.org/MultivariateStats.jl/stable/fa/) or [principal component analysis](https://juliastats.org/MultivariateStats.jl/stable/pca/) model and directly rotate the resulting loading matrix. 

```jldoctest multivariatestats
julia> using FactorRotations, MultivariateStats

julia> data = rand(5, 100);

julia> model = fit(FactorAnalysis, data')
Factor Analysis(indim = 5, outdim = 4)
```
