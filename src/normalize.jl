"""
    kaiser_normalize!(Λ)

Perform an in-place Kaiser normalization of loading matrix `Λ`.

Returns a tuple of a normalized loading matrix and associated weights.

## Examples
```jldoctest
$(DEFINITION_L)
julia> L_norm, weights = kaiser_normalize!(L)
([0.9025386703820677 -0.43060881141120344; 0.8675243270163733 -0.4973947547318815; … ; 0.8009698251886489 0.598704717817778; 0.8891440336983804 0.457627454747389], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

julia> L_norm ≈ L
true
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
julia> L_norm, weights = kaiser_normalize(L)
([0.9025386703820677 -0.43060881141120344; 0.8675243270163733 -0.4973947547318815; … ; 0.8009698251886489 0.598704717817778; 0.8891440336983804 0.457627454747389], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

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
