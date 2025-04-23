"""
    TargetRotation(target::AbstractMatrix; orthogonal = false)

The (partial) target rotation criterion.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The method rotates a factor loading matrix ``Λ ∈ ℝ^{p×k}`` towards the `target` matrix
``T ∈ ℝ^{p×k}``.
For a fully specified `target` matrix (e.g. all its entries are numbers),
the method tries to make the rotated loading matrix as similar to the `target`
as possible:

```math
Q_T(Λ) = \\frac{1}{2} \\left\\| Λ - T \\right\\|^2.
```

When some of the `target` matrix elements are `missing` instead of being finite numbers,
the *partially specified* target rotation is performed.
In that case the criterion becomes

```math
Q_{WT}(Λ) = \\frac{1}{2} \\left\\| (Λ - T) ⋅ W \\right\\|^2,
```

where ``W ∈ ℝ^{p×k}`` is a weight matrix that has *1* for non-missing
elements of *T* and *0* otherwise.

## Examples
### Full target rotation
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> target = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1; 0 1; 0 1];

julia> rotate(L, TargetRotation(target, orthogonal = true))
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.882633  0.258215
 0.922358  0.195806
 0.892467  0.167726
 0.862116  0.233154
 0.252473  0.89669
 0.195508  0.789382
 0.146707  0.726945
 0.260213  0.679549
```

### Partially specified target rotation
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> target = [1 0; missing missing; 1 0; 1 0; 0 1; 0 1; 0 1; 0 1];

julia> rotate(L, TargetRotation(target, orthogonal = true))
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.873299  0.288209
 0.915133  0.227193
 0.886218  0.198109
 0.85365   0.262462
 0.221701  0.90479
 0.168434  0.795599
 0.121793  0.731532
 0.236852  0.68804
```
"""
struct TargetRotation{T,V<:AbstractMatrix,W<:Union{AbstractMatrix, Nothing}} <: RotationMethod{T}
    target::V
    weights::W

    function TargetRotation(target; orthogonal = false)
        T = orthogonal ? Orthogonal : Oblique

        # construct weight matrix assuming 'missing' are unspecified values
        if any(ismissing, target)
            weights = @. !ismissing(target)
            target = coalesce.(target, zero(nonmissingtype(eltype(target))))
        else # no missing values, weights are not required
            weights = nothing
        end

        return new{T,typeof(target),typeof(weights)}(target, weights)
    end
end

function criterion_and_gradient!(∇Q::OptionalGradient, method::TargetRotation, Λ::AbstractMatrix)
    @unpack target, weights = method
    size(target) == size(Λ) ||
        throw(ArgumentError("target matrix and loading matrix must be of equal size"))
    dQ = isnothing(∇Q) ? similar(Λ) : ∇Q
    @. dQ = Λ - target
    if !isnothing(weights)
        @. dQ *= weights
    end
    Q = norm(dQ)^2 / 2
    return Q
end
