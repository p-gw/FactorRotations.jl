"""
    Quartimax
"""
struct Quartimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient(::Quartimax, Λ::AbstractMatrix{<:Real})
    Q = -1 / 4 * sum(x -> x^4, Λ)
    ∇Q = @. -Λ^3
    return Q, ∇Q
end
