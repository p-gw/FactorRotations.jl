"""
    Quartimax()

The Quartimax rotation criterion.

## Details
The Quartimax criterion is a special case of the [`Oblimin`](@ref) rotation criterion with
parameter `gamma = 0`.

## Examples
### Setting up the criterion
```jldoctest
julia> Quartimax()
Quartimax()
```

### Testing equivalence of Quartimax and Oblimin
```jldoctest
$(DEFINITION_L)
julia> L_quartimax = rotate(L, Quartimax());
┌ Info: Rotation algorithm converged after 13 iterations.
│       algorithm: Quartimax
└       criterion: -1.0227347961934468

julia> L_oblimin = rotate(L, Oblimin(gamma = 0, orthogonal = true));
┌ Info: Rotation algorithm converged after 13 iterations.
│       algorithm: Oblimin{Orthogonal, Int64}
└       criterion: 0.1260609090703036

julia> L_quartimax ≈ L_oblimin
true
```

"""
struct Quartimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient(::Quartimax, Λ::AbstractMatrix{<:Real})
    Q = -1 / 4 * sum(x -> x^4, Λ)
    ∇Q = @. -Λ^3
    return Q, ∇Q
end
