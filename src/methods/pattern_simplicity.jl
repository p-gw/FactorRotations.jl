"""
    PatternSimplicity(; orthogonal = false)

The Pattern Simplicity factor rotation criterion.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)
"""
struct PatternSimplicity{RT} <: RotationMethod{RT}
    function PatternSimplicity(; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T}()
    end
end

function criterion(method::PatternSimplicity, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    m = Λsq' * Λsq
    return logdet((diagm(diag(m)))) - logdet(m)
end

function criterion_and_gradient(method::PatternSimplicity, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    m = Λsq' * Λsq
    diag_m = diagm(diag(m))
    Q = logdet(diag_m) - logdet(m)
    ∇Q = -4 * Λ .* (Λsq * (inv(m) - inv(diag_m)))
    return Q, ∇Q
end
