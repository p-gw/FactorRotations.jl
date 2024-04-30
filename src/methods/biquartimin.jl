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

function criterion_and_gradient!(∇Q::OptionalGradient, ::Biquartimin, Λ::AbstractMatrix)
    Λ₂ = @view Λ[:, 2:end]
    Λ₂sq = Λ₂ .^ 2
    part = !isnothing(∇Q) ? @view(∇Q[:, 2:end]) : similar(Λ₂sq)
    part .= sum(Λ₂sq, dims = 2) .- Λ₂sq

    Q = sum(Λ₂sq .* part)

    if !isnothing(∇Q)
        # ∇Q[:, 2:end] === part
        ∇Q[:, 1] .= zero(eltype(∇Q))
        @. part *= 4 * Λ₂
    end
    return Q
end
