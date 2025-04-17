"""
    Varimax()

The *Varimax* rotation criterion.

## Details
The *Varimax* is an orthogonal rotation method that maximizes the column variances of
the loading matrix ``Λ ∈ ℝ^{p × k}``:

```math
Q(Λ) = -\\frac{p}{4} ∑_{j=1}^{k} \\mathrm{Var}_i(Λ_{i, j}).
```

It is a special case of the [`Oblimin`](@ref) rotation criterion with parameter
`gamma = 1`.

## Examples
### Setting up the criterion
```jldoctest
julia> Varimax()
Varimax()
```

### Testing equivalence of Varimax and Oblimin
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> L_varimax = rotate(L, Varimax());

julia> L_oblimin = rotate(L, Oblimin(gamma = 1, orthogonal = true));

julia> loadings(L_varimax) ≈ loadings(L_oblimin)
true
```

"""
struct Varimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, ::Varimax, Λ::AbstractMatrix)
    Λsq = isnothing(∇Q) ? similar(Λ) : ∇Q
    Λsq .= Λ .^ 2
    centercols!(Λsq)
    Q = -norm(Λsq)^2 / 4
    if !isnothing(∇Q)
        @. ∇Q *= -Λ # ∇Q is already centered Λsq
    end
    return Q
end
