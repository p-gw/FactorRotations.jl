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

function criterion_and_gradient(::Varimax, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    colmeans = mean(Λsq, dims = 1)
    Λsq_demeaned = Λsq .- colmeans
    Q = -norm(Λsq_demeaned)^2 / 4
    ∇Q = @. -Λ * Λsq_demeaned
    return Q, ∇Q
end
