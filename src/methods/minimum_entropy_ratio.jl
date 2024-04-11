"""
    MinimumEntropyRatio()

The Minimum Entropy Ratio rotation method.
"""
struct MinimumEntropyRatio <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, ::MinimumEntropyRatio, Λ::AbstractMatrix{T}) where {T}
    p, k = size(Λ)
    Λsq = Λ .^ 2

    total = sum(Λsq)
    colsums = sum(Λsq, dims = 1)

    p₂ = colsums / total

    Q₁ = sum(mxlogx, Λsq / colsums)
    Q₂ = sum(mxlogx, p₂)
    Q = log(Q₁) - log(Q₂)
    isnothing(∇Q) && return Q

    u = Ones(T, p)
    v = Ones(T, k)

    M = u * u'
    R = M * Λsq
    P = Λsq ./ R
    H = @. -(log(P) + 1)
    G₁ = H ./ R - M * (Λsq .* H ./ (R .^ 2))

    h = @. -(log(p₂) + 1)
    α = h * p₂'
    G₂ = u * h / total - α .* u * v'
    @. ∇Q = 2Λ * (G₁ / Q₁ - G₂ / Q₂)

    return Q
end
