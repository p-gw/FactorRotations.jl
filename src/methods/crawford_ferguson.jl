"""
    CrawfordFerguson(; kappa, orthogonal = false)

The family of *Crawford-Ferguson* rotation methods.

## Keyword arguments
- `kappa`: The parameter determining the rotation criterion (see *Details*).
- `orthogonal`: orthogonal: If orthogonal = true an orthogonal rotation is performed, an
                oblique rotation otherwise. (default: `false`)

## Details
The *Crawford-Ferguson* family allows both orthogonal and oblique rotation of the
factor loading matrix with the following criterion:

```math
\\begin{aligned}
Q_{\\mathrm{CF}}(Λ, ϰ) &=
    \\frac{1 - ϰ}{4} \\left⟨Λ², Λ²⋅\\left(1^{k×k} - I\\right) \\right⟩ +
    \\frac{ϰ}{4} \\left⟨Λ², (1^{p×p} - I)⋅Λ²\\right⟩ = \\\\
  &\\hspace{4em} = -\\frac{1 - ϰ}{4} k ∑_{i=1}^p \\mathrm{Var}(Λ²_{i, \\cdot})
    - \\frac{ϰ}{4} p ∑_{j=1}^k \\mathrm{Var}(Λ²_{⋅, j}) = \\\\
  &\\hspace{8em} = \\frac{1 - ϰ}{4} ∑_{i=1}^p \\left(∑_{j=1}^k λ_{i,j}² \\right)²
    + \\frac{ϰ}{4} ∑_{j=1}^k \\left(∑_{i=1}^p λ_{i,j}² \\right)²
    - \\frac{1}{4} ∑_{i,j} λ⁴_{i,j}.
\\end{aligned}
```

If the rotation is *orthogonal*, the *communalities* (``c_i = ∑_{j=1}^k λ_{i,j}²``)
are preserved, and, for a specific `kappa`, Crawford-Ferguson criterion becomes equivalent
to the following rotation methods:
- *ϰ = γ/p*: [`Oblimin(gamma = γ, orthogonal = true)`](@ref Oblimin)
- *ϰ = 0*: [`Quartimax`](@ref)
- *ϰ = 1/p*: [`Varimax`](@ref)
- *ϰ = k/2p*: [`Equamax`](@ref)
- *ϰ = (k - 1)/(p + k - 2)*: [`Parsimax`](@ref)
- *ϰ = 1*: Factor parsimony

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
