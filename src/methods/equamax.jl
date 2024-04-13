"""
    Equamax()

*Equamax* is an orthogonal rotation method, which is equivalent to [`Oblimin`](@ref)
rotation of *p × k* loadings matrix *Λ* with ``\\gamma = \\frac{p}{2}``.
"""
struct Equamax <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, ::Equamax, Λ::AbstractMatrix{<:Real})
    p, k = size(Λ)
    weighted_sums_criterion_and_gradient!(∇Q, Λ, 1 - k / (2p), k / (2p))
end
