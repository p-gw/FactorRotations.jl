"""
    CrawfordFerguson(; kappa, orthogonal = false)

The family of Crawford-Ferguson rotation methods.

## Keyword arguments
- `kappa`: The parameter determining the rotation criterion (see Details).
- `orthogonal`: orthogonal: If orthogonal = true an orthogonal rotation is performed, an
                oblique rotation otherwise. (default: `false`)

## Details
The Crawford-Ferguson family allows both orthogonal and oblique rotation of the
`p`-by`k`-factor loading matrix. If orthogonal rotation is performed, Crawford-Ferguson is
equivalent to Oblimin rotation given the following values for `kappa`:

- `kappa = 0` → [`Quartimax`](@ref)
- `kappa = 1/p` → [`Varimax`](@ref)
- `kappa = k/2p` → Equamax
- `kappa = (k - 1)/(p + k - 2)` → Parsimax
- `kappa = 1` → Factor parsimony

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

function criterion_and_gradient!(∇Q::OptionalGradient, method::CrawfordFerguson, Λ::AbstractMatrix{T}) where {T}
    @unpack κ = method
    p, k = size(Λ)

    N = Ones(k, k) - I(k)
    M = Ones(p, p) - I(p)

    Λsq = Λ .^ 2

    ΛsqN = Λsq * N
    MΛsq = M * Λsq

    Q = ((1 - κ)/4) * tr(Λsq' * ΛsqN) + (κ/4) * tr(Λsq' * MΛsq)
    if !isnothing(∇Q)
        @. ∇Q = (1 - κ) * Λ * ΛsqN + κ * Λ * MΛsq
    end
    return Q
end
