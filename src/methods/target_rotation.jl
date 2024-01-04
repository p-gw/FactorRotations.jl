"""
    TargetRotation
"""
struct TargetRotation{T,V<:AbstractMatrix} <: RotationMethod{T}
    H::V
    W::BitMatrix
    function TargetRotation(target; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique

        # construct weight matrix assuming 'missing' and/or 'nothing' are unspecified
        W = @. !ismissing(target)
        H = coalesce.(target, 0)

        return new{T,typeof(target)}(H, W)
    end
end

function criterion_and_gradient(method::TargetRotation, Λ::AbstractMatrix)
    @unpack H, W = method
    size(H) == size(Λ) ||
        throw(ArgumentError("target matrix and loading matrix must be of equal size"))
    ∇Q = @. W * (Λ - H)
    Q = norm(∇Q)^2 / 2
    return Q, ∇Q
end
