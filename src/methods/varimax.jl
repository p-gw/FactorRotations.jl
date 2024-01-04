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
```jldoctest
$(DEFINITION_L)
julia> L_varimax = rotate(L, Varimax());
┌ Info: Rotation algorithm converged after 9 iterations.
│       algorithm: Varimax
└       criterion: -0.4515671564134383

julia> L_oblimin = rotate(L, Oblimin(gamma = 1, orthogonal = true));
┌ Info: Rotation algorithm converged after 9 iterations.
│       algorithm: Oblimin{Orthogonal, Int64}
└       criterion: -0.4149267008747196

julia> L_varimax ≈ L_oblimin
true
```

"""
struct Varimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient(::Varimax, Λ::AbstractMatrix{T}) where {T<:Real}
    Λsq = Λ .^ 2
    colmeans = mean(Λsq, dims = 1)
    Λsq_demeaned = Λsq .- colmeans
    Q = -norm(Λsq_demeaned)^2 / 4
    ∇Q = @. -Λ * Λsq_demeaned
    return Q, ∇Q
end
