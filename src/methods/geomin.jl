"""
    Geomin(epsilon = 0.01)

The Geomin rotation method.

## Keyword arguments
- `epsilon`: A small constant to deal with zero loadings.
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
    Λsq = Λ .^ 2 .+ ε
    p, k = size(Λ)
    u = Ones(T, p)
    v = Ones(T, k)

    part = exp.((log.(Λsq) * v) ./ k)
    Q = u' * part
    if !isnothing(∇Q)
        mul!(∇Q, part, v')
        ∇Q .*= (2 / k) .* inv.(Λsq) .* Λ
    end

    return Q
end
