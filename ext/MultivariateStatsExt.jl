module MultivariateStatsExt

using FactorRotations, MultivariateStats

"""
    rotate(model, method; kwargs...)

Perform a rotation of the models loading matrix using rotation `method`, where model is a
fitted `MultivariateStats.FactorAnalysis`, `MultivariateStats.PCA`, or
`MultivariateStats.PPCA` model.

For a list of available keyword arguments see [`rotate`](@ref).
"""
function FactorRotations.rotate(model::Union{PCA,PPCA,FactorAnalysis}, method; kwargs...)
    L = loadings(model)
    return rotate(L, method; kwargs...)
end

"""
    rotate!(model, method; kwargs...)

Perform an in-place rotation of the models loading matrix using rotation `method`, where model is a
fitted `MultivariateStats.FactorAnalysis`, `MultivariateStats.PCA`, or
`MultivariateStats.PPCA` model.

For a list of available keyword arguments see [`rotate`](@ref).
"""
function FactorRotations.rotate!(model::Union{PCA,PPCA,FactorAnalysis}, method; kwargs...)
    L = loadings(model)
    return rotate!(L, method; kwargs...)
end

end
