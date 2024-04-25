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
    p, k = size(Λ)
    Λsq = Λ .^ 2

    total = sum(Λsq)
    Λsq ./= total
    rowsums = sum(Λsq, dims = 2)
    colsums = sum(Λsq, dims = 1)

    Q = -log(k) + sum(mxlogx, Λsq) - sum(mxlogx, rowsums) -
        sum(mxlogx, colsums)
    isnothing(∇Q) && return Q

    u = Ones(T, p)
    v = Ones(T, k)

    H = @. -(log(Λsq) + 1)
    α = u' * (Λsq .* H) * v / total
    G₀ = H / total - α * u * v'

    h₁ = -(log.(Λsq * v) .+ 1)
    α₁ = v' * Λsq' * h₁ / total
    G₁ = h₁ * v' / total - α₁ * u * v'

    h₂ = -(log.(u' * Λsq) .+ 1)
    α₂ = h₂ * Λsq' * u / total
    G₂ = u * h₂ / total - α₂ .* u * v'

    @. ∇Q = 2Λ * (G₀ - G₁ - G₂)
    return Q
end
