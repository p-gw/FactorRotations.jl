"""
minetropy
"""
struct MinimumEntropy <: RotationMethod{Orthogonal} end

function criterion_and_gradient(::MinimumEntropy, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    logΛsq = log.(Λsq)
    Q = -tr(Λsq' * logΛsq) / 2
    ∇Q = -Λ .* logΛsq - Λ
    return Q, ∇Q
end
