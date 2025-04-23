"""
    Parsimax()

*Parsimax* is an orthogonal rotation method, which is equivalent to [`CrawfordFerguson`](@ref)
rotation of *p × k* loadings matrix *Λ* with ``ϰ = \\frac{k - 1}{p + k - 2}``.
"""
struct Parsimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, ::Parsimax, Λ::AbstractMatrix{<:Real})
    p, k = size(Λ)
    weighted_sums_criterion_and_gradient!(∇Q, Λ, (p - 1) / (p + k - 2), (k - 1) / (p + k - 2))
end
