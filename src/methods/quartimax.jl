"""
    Quartimax()

The *Quartimax* rotation criterion.

## Details
The *Quartimax* criterion is a special case of the [`Oblimin`](@ref) rotation criterion with
parameter `gamma = 0`.

## Examples
### Setting up the criterion
```jldoctest
julia> Quartimax()
Quartimax()
```

### Testing equivalence of Quartimax and Oblimin
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> L_quartimax = rotate(L, Quartimax());

julia> L_oblimin = rotate(L, Oblimin(gamma = 0, orthogonal = true));

julia> loadings(L_quartimax) ≈ loadings(L_oblimin)
true
```

"""
struct Quartimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient!(∇Q::OptionalGradient, method::Quartimax, Λ::AbstractMatrix)
    if !isnothing(∇Q)
        @. ∇Q = -Λ^3
    end
    return -sum(x -> x^4, Λ) / 4
end
