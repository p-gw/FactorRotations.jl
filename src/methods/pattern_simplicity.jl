"""
    PatternSimplicity(; orthogonal = false)

The *Pattern Simplicity* factor rotation criterion.

## Details

```math
Q_{\\mathrm{PS}}(Λ) = \\log \\left\\| \\mathrm{diag}\\left(\\left(Λ²\\right)ᵀ ⋅ Λ²\\right) \\right\\| -
                      \\log \\left\\| \\left(Λ²\\right)ᵀ ⋅ Λ² \\right\\|,
```

where ``Λ²`` is the matrix of squared loadings.

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

function criterion_and_gradient!(∇Q::OptionalGradient, method::PatternSimplicity, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    m = Λsq' * Λsq
    diag_m = diagm(diag(m))
    if !isnothing(∇Q)
        ∇Q .= Λ .* (Λsq * (inv(m) - inv(diag_m)))
        lmul!(-4, ∇Q)
    end
    return logdet(diag_m) - logdet(m)
end
