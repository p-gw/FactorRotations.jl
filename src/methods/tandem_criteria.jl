"""
    TandemCriterionI()

The first criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionI <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q, method::TandemCriterionI, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    part = ΛxΛ .^ 2 * Λsq
    if !isnothing(∇Q)
        ∇Q .= Λ .* part .+ (ΛxΛ .* (Λsq * Λsq')) * Λ
        lmul!(-4, ∇Q)
    end
    return -tr(Λsq' * part)
end

"""
    TandemCriterionI()

The second criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionII <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q, method::TandemCriterionII, Λ::AbstractMatrix)
    p, k = size(Λ)
    u = Ones(p)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    ΛxΛsq = ΛxΛ .^ 2
    part = (u * u' - ΛxΛsq) * Λsq
    if !isnothing(∇Q)
        ∇Q .= Λ .* part
        mul!(∇Q, ΛxΛ .* (Λsq * Λsq'), Λ, -4, 4)
    end
    return tr(Λsq' * part)
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
