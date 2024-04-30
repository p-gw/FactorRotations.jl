"""
    Oblimin(; gamma, orthogonal = false)

The family of *Oblimin* rotation methods.

## Keyword arguments
- `gamma`: The shape parameter determining the rotation criterion (see *Details*).
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The *Oblimin* rotation family allows orthogonal as well as oblique rotation of the *p*×*k* factor
loading matrix. If orthogonal rotation is performed, *Oblimin* is equivalent to the following
rotation methods given a value for `gamma`:
- *γ = p×κ*: [`CrawfordFerguson(kappa = κ, orthogonal = true)`](@ref CrawfordFerguson)
- *γ = 0*: [`Quartimax`](@ref)
- *γ = 1/2*: [`Biquartimax`](@ref)
- *γ = 1*: [`Varimax`](@ref)
- *γ = k/2*: [`Equamax`](@ref)
- *γ = p×(k - 1)/(p + k - 2)*: [`Parsimax`](@ref)

For oblique rotation *Oblimin* is equivalent to the following rotation methods:

- *γ = 0*: *Quartimin*
- *γ = 1/2*: [`Biquartimin`](@ref)

## Examples
```jldoctest
julia> Oblimin(gamma = 0.5)
Oblimin{Oblique, Float64}(0.5)

julia> Oblimin(gamma = 1, orthogonal = true)
Oblimin{Orthogonal, Int64}(1)
```
"""
struct Oblimin{T,V} <: RotationMethod{T}
    γ::V
    function Oblimin(; gamma, orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T,typeof(gamma)}(gamma)
    end
end

criterion_and_gradient!(∇Q::OptionalGradient, method::Oblimin, Λ::AbstractMatrix{<:Real}) =
    weighted_sums_criterion_and_gradient!(∇Q, Λ, 1 - method.γ, method.γ / size(Λ, 1))
