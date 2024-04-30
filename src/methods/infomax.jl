"""
    Infomax(; orthogonal = false)

The Infomax rotation method.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)
"""
struct Infomax{T} <: RotationMethod{T}
    function Infomax(; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T}()
    end
end

function criterion_and_gradient!(∇Q::OptionalGradient, ::Infomax, Λ::AbstractMatrix{T}) where {T}
    k = size(Λ, 2)
    Λsq = Λ .^ 2

    total = sum(Λsq)
    Λsq ./= total
    rowsums = sum(Λsq, dims = 2)
    colsums = sum(Λsq, dims = 1)

    Q = -log(k) + sum(mxlogx, Λsq) - sum(mxlogx, rowsums) -
        sum(mxlogx, colsums)
    isnothing(∇Q) && return Q

    H = @. -(log(Λsq) + 1)
    G₀ = @. (H - dot(Λsq, H))

    h₁ = @. -(log(rowsums) + 1)
    G₁ = @. (h₁ - dot(rowsums, h₁))

    h₂ = @. -(log(colsums) + 1)
    G₂ = @. (h₂ - dot(h₂, colsums))

    @. ∇Q = (2/total) * Λ * (G₀ - G₁ - G₂)
    return Q
end
