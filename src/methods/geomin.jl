"""
    Geomin(epsilon = 0.01)

The *Geomin* rotation method.

## Details
The *Geomin* is an oblique rotation method that minimizes the
sum of column-wise geometric means of the squared loadings:

```math
Q_{\\mathrm{Geomin}}(Λ, ε) =
    ∑_{i = 1}^p \\left(\\prod_{j = 1}^k λ²_{i, j} + ε \\right)^{\\frac{1}{k}}.
```

## Keyword arguments
- `epsilon`: A small non-negative constant to deal with zero loadings.
"""
struct Geomin{T} <: RotationMethod{Oblique}
    ε::T
    function Geomin(; epsilon::T = 0.01) where {T<:Real}
        epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
        return new{T}(epsilon)
    end
end

function criterion_and_gradient!(∇Q::OptionalGradient, method::Geomin, Λ::AbstractMatrix{T}) where {T}
    @unpack ε = method
    Λsq = !isnothing(∇Q) ? ∇Q : similar(Λ)
    Λsq .= Λ .^ 2 .+ ε
    k = size(Λ, 2)

    part = exp.(sum(log.(Λsq), dims=2) ./ k)
    Q = sum(part)
    if !isnothing(∇Q)
        # ∇Q === Λsq
        ∇Q .= (2 / k) .* Λ ./ Λsq .* part
    end

    return Q
end
