# Working with MultivariateStats.jl

FactorRotations.jl provides direct support for models fitted by [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl). Specifically, you can fit a [factor analysis](https://juliastats.org/MultivariateStats.jl/stable/fa/) or [principal component analysis](https://juliastats.org/MultivariateStats.jl/stable/pca/) model and directly rotate the resulting loading matrix. 

To user FactorRotations.jl with MultivariateStats.jl we first need to fit a factor analysis or principal component analysis model.

```jldoctest multivariatestats
julia> using FactorRotations, MultivariateStats

julia> data = rand(5, 100);  # TODO: use some real data set

julia> model = fit(FactorAnalysis, data')
Factor Analysis(indim = 5, outdim = 4)
```

Analogous to rotating raw loading matrices, there are two ways to rotate a MultivariateStats.jl solution: regular and in-place. 
For both we can pass our `model` to [`rotate`](@ref) or [`rotate!`](@ref) respectively.

```jldoctest multivariatestats
julia> r = rotate(model, Geomin())
FactorRotation{Float64} with loading matrix:
...
```

Similarly we can use [`rotate!`](@ref) to change the loading matrix of `model`.

```jldoctest multivariatestats
julia> rotate!(model, Geomin())
...
```
