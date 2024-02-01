# Working with MultivariateStats.jl

FactorRotations.jl provides direct support for models fitted by [MultivariateStats.jl](https://github.com/JuliaStats/MultivariateStats.jl). Specifically, you can fit a [factor analysis](https://juliastats.org/MultivariateStats.jl/stable/fa/) or [principal component analysis](https://juliastats.org/MultivariateStats.jl/stable/pca/) model and directly rotate the resulting loading matrix. 

To user FactorRotations.jl with MultivariateStats.jl we first need to fit a factor analysis or principal component analysis model.
In this tutorial we will be using the `bfi` data provided by the [`psych`](https://cran.r-project.org/web/packages/psych/psych.pdf) R package.
It contains 25 self report items concerning personality. 

!!! note
    For performance reasons we just use the first 200 observations of the dataset.

!!! warning
    Missing values are dropped for purposes of this tutorial. 
    Properly handle missing data in your own analysis!

```jldoctest multivariatestats
julia> using FactorRotations, MultivariateStats, RDatasets

julia> data = dataset("psych", "bfi")[:, 2:26] |> dropmissing!;

julia> data = Matrix(data[1:200, :]);  # use just the first 200 observations

julia> model = fit(FactorAnalysis, data', maxoutdim = 5)
Factor Analysis(indim = 25, outdim = 5)
```

After fitting the model the loadings could be extracted using `MultivariateStats.loadings` and then rotated,

```jldoctest multivariatestats
julia> raw_loadings = MultivariateStats.loadings(model);

julia> rotated_loadings = FactorRotations.loadings(rotate(raw_loadings, Geomin()));

```

However, FactorAnalysis.jl provides convenience functions to pass `model` directly.
Analogous to rotating raw loading matrices, there are two ways to rotate a MultivariateStats.jl solution: regular and in-place. 
For both we can pass our `model` to [`rotate`](@ref) or [`rotate!`](@ref) respectively.

Using [`rotate`](@ref),

```jldoctest multivariatestats
julia> rot = rotate(model, Geomin());

julia> FactorRotations.loadings(rot) == rotated_loadings
true
```

Similarly we can use [`rotate!`](@ref) to change the loading matrix of `model` in-place.

```jldoctest multivariatestats
julia> rotate!(model, Geomin());

julia> MultivariateStats.loadings(model) == rotated_loadings
true
```
