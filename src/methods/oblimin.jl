"""
    Oblimin(; gamma, orthogonal = false)

The family of Oblimin rotation methods.

## Keyword arguments
- `gamma`: The shape parameter determining the rotation criterion (see Details).
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The Oblimin rotation family allow orthogonal as well as oblique rotation of the factor
loading matrix. If orthogonal rotation is performed, Oblimin is equivalent to the following
rotation methods given a value for `gamma`:

- `gamma = 0` → [`Quartimax`](@ref)
- `gamma = 0.5` → [`Biquartimax`](@ref)
- `gamma = 1` → [`Varimax`](@ref)
- `gamma = p/2` → Equamax

For oblique rotation Oblimin is equivalent to the following rotation methods:

- `gamma = 0` → Quartimin
- `gamma = 0.5` → Biquartimin

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
