"""
    MinimumEntropy()

The *Minimum Entropy* rotation method.

## Details

The *Minimum Entropy* rotation method minimizes the *entropy* of the squared loadings:

```math
Q_{\\mathrm{MinEnt}}(Λ) = -\\frac{1}{2} ∑_{i,j} λ_{i,j}² \\log λ_{i,j}².
```

## See also

[`MinimumEntropyRatio`](@ref)
"""
struct MinimumEntropy <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, ::MinimumEntropy, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    mlogΛsq = !isnothing(∇Q) ? ∇Q : similar(Λ)
    @. mlogΛsq = -log(Λsq)
    Q = dot(Λsq, mlogΛsq) / 2
    if !isnothing(∇Q)
        # ∇Q === mlogΛsq
        ∇Q .-= one(eltype(∇Q))
        ∇Q .*= Λ
    end
    return Q
end
