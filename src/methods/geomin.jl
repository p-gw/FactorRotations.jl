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

function criterion(method::Geomin, Λ::AbstractMatrix{T}) where {T}
    @unpack ε = method
    Λsq = Λ .^ 2 .+ ε
    p, k = size(Λ)
    u = Ones(T, p)
    v = Ones(T, k)
    Q = u' * exp.((log.(Λsq) ./ k * v))
    return Q
end

function criterion_and_gradient(method::Geomin, Λ::AbstractMatrix{T}) where {T}
    @unpack ε = method
    Λsq = Λ .^ 2 .+ ε
    p, k = size(Λ)
    u = Ones(T, p)
    v = Ones(T, k)

    part = exp.((log.(Λsq) ./ k * v))
    Q = u' * part
    ∇Q = (2 ./ k) .* (1 ./ Λsq) .* Λ .* (part * v')

    return Q, ∇Q
end
