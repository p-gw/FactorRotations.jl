"""
    Biquartimin

The Biquartimin rotation criterion.

# Keyword arguments
- `orthogonal`: orthogonal: If orthogonal = true an orthogonal rotation is performed, an
                oblique rotation otherwise. (default: `false`)
"""
struct Biquartimin{RT} <: RotationMethod{RT}
    function Biquartimin(; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T}()
    end
end

function criterion_and_gradient(::Biquartimin, Λ::AbstractMatrix)
    p, k = size(Λ)
    n = k - 1

    Λ₂ = @view Λ[:, 2:end]
    Λ₂sq = Λ₂ .^ 2
    N = Ones(n, n) - I(n)

    Q = sum(Λ₂sq .* (Λ₂sq * N))

    ∇Q = zeros(size(Λ))
    ∇Q[:, 2:end] .= 4 * Λ₂ .* (Λ₂sq * N)
    return Q, ∇Q
end
