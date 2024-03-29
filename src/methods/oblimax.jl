"""
    Oblimax(; orthogonal = false)

The Oblimax rotation method.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The Oblimax rotation method is equivalent to [`Quartimax`](@ref) for orthogonal rotation.

## Examples
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> L_oblimax = rotate(L, Oblimax(orthogonal = true));

julia> L_quartimax = rotate(L, Quartimax());

julia> isapprox(loadings(L_oblimax), loadings(L_quartimax), atol = 1e-6)
true
```
"""
struct Oblimax{T} <: RotationMethod{T}
    function Oblimax(; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique
        return new{T}()
    end
end

function criterion_and_gradient(::Oblimax, Λ::AbstractMatrix)
    sqnorm_Λsq = norm(Λ .^ 2) .^ 2
    norm_Λ = norm(Λ)

    K = sqnorm_Λsq / norm_Λ^4
    Q = -log(K)
    ∇Q = -4Λ .^ 3 / sqnorm_Λsq + 4Λ / norm_Λ^2
    return Q, ∇Q
end
