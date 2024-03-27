"""
    kaiser_normalize!(Λ)

Perform an in-place Kaiser normalization of loading matrix `Λ`.

Returns a tuple of a normalized loading matrix and associated weights.

## Examples
```jldoctest
$(DEFINITION_L)
julia> L_norm, weights = kaiser_normalize!(L);

julia> L_norm
8×2 Matrix{Float64}:
 0.902539  -0.430609
 0.867524  -0.497395
 0.855641  -0.517569
 0.89353   -0.449004
 0.84375    0.536737
 0.826331   0.563184
 0.80097    0.598705
 0.889144   0.457627

julia> weights
8-element Vector{Float64}:
 0.9196281857359527
 0.9429130394686458
 0.9080908544853868
 0.8930873417533137
 0.9315556880831118
 0.8132330539273475
 0.7416009708731509
 0.7276661322337326

julia> L_norm ≈ L
true

```
"""
function kaiser_normalize!(Λ::AbstractMatrix)
    weights = communalities(Λ)
    Λ ./= weights
    return Λ, weights
end

"""
    kaiser_normalize(Λ)

Perform a Kaiser normalization of loading matrix `Λ`.

Returns a tuple of a normalized loading matrix and associated weights.

## Examples
```jldoctest
$(DEFINITION_L)
julia> L_norm, weights = kaiser_normalize(L);

julia> L_norm
8×2 Matrix{Float64}:
 0.902539  -0.430609
 0.867524  -0.497395
 0.855641  -0.517569
 0.89353   -0.449004
 0.84375    0.536737
 0.826331   0.563184
 0.80097    0.598705
 0.889144   0.457627

julia> weights
8-element Vector{Float64}:
 0.9196281857359527
 0.9429130394686458
 0.9080908544853868
 0.8930873417533137
 0.9315556880831118
 0.8132330539273475
 0.7416009708731509
 0.7276661322337326

```
"""
kaiser_normalize(Λ) = kaiser_normalize!(copy(Λ))

"""
    kaiser_denormalize!(Λ, weights)

Undo a Kaiser normalization of normalized `Λ` in-place given `weights`.

## Examples
```
$(DEFINITION_L)
julia> L_orig = copy(L);

julia> _, weights = kaiser_normalize!(L);

julia> kaiser_denormalize!(L, weights)
8×2 Matrix{Float64}:
 0.83   -0.396
 0.818  -0.469
 0.777  -0.47
 0.798  -0.401
 0.786   0.5
 0.672   0.458
 0.594   0.444
 0.647   0.333

julia> L ≈ L_orig
true

```
"""
function kaiser_denormalize!(Λ::AbstractMatrix, weights::AbstractVector)
    Λ .*= weights
    return Λ
end

"""
    kaiser_denormalize(Λ, weights)

Undo a Kaiser normalization of normalized loading matrix `Λ` given `weights`.

## Examples
```jldoctest
$(DEFINITION_L)
julia> L_norm, weights = kaiser_normalize(L);

julia> L_denorm = kaiser_denormalize(L_norm, weights)
8×2 Matrix{Float64}:
 0.83   -0.396
 0.818  -0.469
 0.777  -0.47
 0.798  -0.401
 0.786   0.5
 0.672   0.458
 0.594   0.444
 0.647   0.333

julia> L ≈ L_denorm
true

```
"""
kaiser_denormalize(Λ, weights) = kaiser_denormalize!(copy(Λ), weights)

communalities(Λ) = norm.(eachrow(Λ))
