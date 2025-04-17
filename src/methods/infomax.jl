"""
    Infomax(; orthogonal = false)

The *Infomax* rotation method.

## Details

The *Infomax* method maximizes the mutual information between the variables and factors:

```math
\\begin{aligned}
Q_{\\mathrm{Infomax}}(Λ) =& \\left(∑_{i,j} λ_{i,j}²\\right)^{-1} \\left(
    - ∑_{i,j} λ_{i,j}² \\log λ_{i,j}²
    + ∑_{i=1}^p \\left(∑_{j=1}^k λ_{i,j}²\\right) \\log ∑_{j=1}^k λ_{i,j}² + \\right. \\\\
    & \\hspace{8em} \\left.
    + ∑_{j=1}^k \\left(∑_{i=1}^p λ_{i,j}²\\right) \\log ∑_{i=1}^p λ_{i,j}²
    - \\log k
    \\right).
\\end{aligned}
```

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
