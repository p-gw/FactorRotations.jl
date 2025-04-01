"""
    Biquartimax()

The *Biquartimax* rotation method.

## Details
The *Biquartimax* rotation method is a special case of the [`Oblimin`](@ref) rotation with
parameters `gamma = 0.5` and `orthogonal = true`.

## Examples
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> L_biquartimax = rotate(L, Biquartimax());

julia> L_oblimin = rotate(L, Oblimin(gamma = 0.5, orthogonal = true));

julia> loadings(L_biquartimax) ≈ loadings(L_oblimin)
true
```
"""
Biquartimax() = Oblimin(gamma = 0.5, orthogonal = true)
