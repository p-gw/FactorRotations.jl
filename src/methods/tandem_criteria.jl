"""
    TandemCriterionI()

The first criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionI <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, method::TandemCriterionI, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    part = !isnothing(∇Q) ? ∇Q : similar(Λsq)
    mul!(part, ΛxΛ.^2, Λsq)
    Q = -dot(Λsq, part)
    if !isnothing(∇Q)
        ∇Q .*= Λ # ∇Q === part
        mul!(∇Q, ΛxΛ .* (Λsq * Λsq'), Λ, -4, -4)
    end
    return Q
end

"""
    TandemCriterionI()

The second criterion of the tandem criteria factor rotation method.
"""
struct TandemCriterionII <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, method::TandemCriterionII, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    ΛxΛ = Λ * Λ'
    part = !isnothing(∇Q) ? ∇Q : similar(Λsq)
    mul!(part, (1 .- ΛxΛ.^2), Λsq)
    Q = dot(Λsq, part)
    if !isnothing(∇Q)
        ∇Q .*= Λ # ∇Q === part
        mul!(∇Q, ΛxΛ .* (Λsq * Λsq'), Λ, -4, 4)
    end
    return Q
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
