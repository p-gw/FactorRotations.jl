"""
    Varimax()

The Varimax rotation criterion.

## Details
The Varimax criterion is a special case of the Oblimin rotation criterion with parameter
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

function criterion_and_gradient!(∇Q, ::Varimax, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    centercols!(Λsq)
    if !isnothing(∇Q)
        @. ∇Q = -Λ * Λsq
    end
    return -norm(Λsq)^2 / 4
end
