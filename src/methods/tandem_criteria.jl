"""
    TandemCriterionI()

The first criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionI <: RotationMethod{Orthogonal} end

function criterion(method::TandemCriterionI, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    return -tr(Λsq' * ((Λ * Λ') .^ 2 * Λsq))
end

function criterion_and_gradient(method::TandemCriterionI, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    part = ΛxΛ .^ 2 * Λsq
    Q = -tr(Λsq' * part)
    ∇Q = -4 * Λ .* part - 4 * (ΛxΛ .* (Λsq * Λsq')) * Λ
    return Q, ∇Q
end

"""
    TandemCriterionI()

The second criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionII <: RotationMethod{Orthogonal} end

function criterion(method::TandemCriterionII, Λ::AbstractMatrix)
    p, k = size(Λ)
    u = Ones(p)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    ΛxΛsq = ΛxΛ .^ 2
    part = (u * u' - ΛxΛsq) * Λsq
    Q = tr(Λsq' * part)
    return Q
end

function criterion_and_gradient(method::TandemCriterionII, Λ::AbstractMatrix)
    p, k = size(Λ)
    u = Ones(p)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    ΛxΛsq = ΛxΛ .^ 2
    part = (u * u' - ΛxΛsq) * Λsq
    Q = tr(Λsq' * part)

    ∇Q = 4 * Λ .* part - 4 * (ΛxΛ .* (Λsq * Λsq')) * Λ
    return Q, ∇Q
end

"""
    TandemCriteria(; keep)

The tandem criteria rotation method.

## Keyword arguments
- `keep`: The number of factors to keep for the second tandem criterion.
"""
struct TandemCriteria <: RotationMethod{Orthogonal}
    keep::Int
    function TandemCriteria(; keep)
        keep > 1 || throw(ArgumentError("must keep more than 1 factor for rotation"))
        return new(keep)
    end
end
