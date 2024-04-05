"""
    MinimumEntropy()

The Minimum Entropy rotation method.
"""
struct MinimumEntropy <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q, ::MinimumEntropy, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    logΛsq = log.(Λsq)
    Q = -tr(Λsq' * logΛsq) / 2
    if !isnothing(∇Q)
        @. ∇Q = (-1 - logΛsq) * Λ
    end
    return Q
end
