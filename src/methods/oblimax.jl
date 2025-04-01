"""
    Oblimax(; orthogonal = false)

The *Oblimax* rotation method.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The *Oblimax* rotation method is equivalent to [`Quartimax`](@ref) for orthogonal rotation.

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

function criterion_and_gradient!(∇Q::OptionalGradient, ::Oblimax, Λ::AbstractMatrix)
    sqnorm_Λsq = sum(x -> x^4, Λ)
    sqnorm_Λ = norm(Λ)^2

    K = sqnorm_Λsq / sqnorm_Λ^2
    Q = -log(K)
    if !isnothing(∇Q)
        ∇Q .= Λ .^ 3
        axpby!(4/sqnorm_Λ, Λ, -4/sqnorm_Λsq, ∇Q)
    end
    return Q
end
