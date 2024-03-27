"""
    reflect!(r::FactorRotation)

Modify `r` in-place by swapping signs of the loading matrix `r.L` such that the sum of each
column is positive. The rotation matrix `r.T` and factor correlation matrix `r.phi` are
updated accordingly.

## Examples
```jldoctest
$(DEFINITION_L)

julia> r = rotate(L, Varimax());

julia> r_reflected = reflect!(r)
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.886061  0.246196
 0.924934  0.183253
 0.894664  0.155581
 0.865205  0.221416
 0.264636  0.893176
 0.206218  0.786653
 0.156572  0.724884
 0.269424  0.67595

julia> r == r_reflected
true
```
"""
function reflect!(r::FactorRotation)
    @unpack L, T, phi = r
    v = reflect_cols(L)
    L .*= v
    T .*= v
    phi .= T' * T
    return r
end

"""
    reflect(r::FactorRotation)

Return a new [`FactorRotation`](@ref) with a modified loading matrix such that the sum of
each column is positive. The rotation matrix and factor correlation matrix are updated
accordingly.

## Examples
```jldoctest
$(DEFINITION_L)

julia> r = rotate(L, Varimax());

julia> reflect(r)
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.886061  0.246196
 0.924934  0.183253
 0.894664  0.155581
 0.865205  0.221416
 0.264636  0.893176
 0.206218  0.786653
 0.156572  0.724884
 0.269424  0.67595

```
"""
reflect(r::FactorRotation) = reflect!(deepcopy(r))

function reflect_cols(m::AbstractMatrix)
    colsums = sum(m, dims = 1)
    return @. ifelse(colsums < 0, -1, 1)
end
