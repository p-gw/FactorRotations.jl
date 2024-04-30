"""
    CrawfordFerguson(; kappa, orthogonal = false)

The family of Crawford-Ferguson rotation methods.

## Keyword arguments
- `kappa`: The parameter determining the rotation criterion (see *Details*).
- `orthogonal`: orthogonal: If orthogonal = true an orthogonal rotation is performed, an
                oblique rotation otherwise. (default: `false`)

## Details
The Crawford-Ferguson family allows both orthogonal and oblique rotation of the
*p*×*k* factor loading matrix. If orthogonal rotation is performed, Crawford-Ferguson with
a specific value for `kappa` is equivalent to the following rotation methods:
- *κ = γ/p*: [`Oblimin(gamma = γ, orthogonal = true)`](@ref Oblimin)
- *κ = 0*: [`Quartimax`](@ref)
- *κ = 1/p*: [`Varimax`](@ref)
- *κ = k/2p*: [`Equamax`](@ref)
- *κ = (k - 1)/(p + k - 2)*: [`Parsimax`](@ref)
- *κ = 1*: Factor parsimony

## Examples
```jldoctest
julia> CrawfordFerguson(kappa = 0, orthogonal = true)
CrawfordFerguson{Orthogonal, Int64}(0)

julia> CrawfordFerguson(kappa = 1/2, orthogonal = false)
CrawfordFerguson{Oblique, Float64}(0.5)
```
"""
struct CrawfordFerguson{T,V} <: RotationMethod{T}
    κ::V
    function CrawfordFerguson(; kappa, orthogonal = false)
        0 <= kappa <= 1 || throw(ArgumentError("kappa must be between 0 and 1"))
        T = orthogonal ? Orthogonal : Oblique
        return new{T,typeof(kappa)}(kappa)
    end
end

criterion_and_gradient!(∇Q::OptionalGradient, method::CrawfordFerguson, Λ::AbstractMatrix{<:Real}) =
    weighted_sums_criterion_and_gradient!(∇Q, Λ, 1 - method.κ, method.κ)
