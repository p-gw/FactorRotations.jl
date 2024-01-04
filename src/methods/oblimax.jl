"""
    Oblimax
"""
struct Oblimax{T} <: RotationMethod{T}
    function Oblimax(; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T}()
    end
end

function criterion_and_gradient(::Oblimax, Λ::AbstractMatrix{T}) where {T<:Real}
    sqnorm_Λsq = norm(Λ .^ 2) .^ 2
    norm_Λ = norm(Λ)

    K = sqnorm_Λsq / norm_Λ^4
    Q = -log(K)
    ∇Q = -4Λ .^ 3 / sqnorm_Λsq + 4Λ / norm_Λ^2
    return Q, ∇Q
end
