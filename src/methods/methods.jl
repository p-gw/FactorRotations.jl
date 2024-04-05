"""
    RotationMethod{T<:RotationType}

An abstract type representing a factor rotation method.

Each implementation of `M <: RotationMethod` must implement at least one of the following methods:
- [`criterion`](@ref)
- [`criterion_and_gradient`](@ref)

If only [`criterion`](@ref) is implemented, gradients are calculated by automatic differentiation.
"""
abstract type RotationMethod{T<:RotationType} end

"""
    criterion(method::RotationMethod, Λ::Abstractmatrix{<:Real})

Calculate the criterion of a given `method` with respect to the factor loading matrix `Λ`.
"""
criterion(method::RotationMethod, Λ::AbstractMatrix{<:Real}) = criterion_and_gradient!(nothing, method, Λ)

"""
    criterion_and_gradient!(∇Q::Union{AbstractMatrix{<:Real}, Nothing},
                            method::RotationMethod, Λ::AbstractMatrix{<:Real})

Calculate the quality criterion *Q* and the gradient of a given `method`
with respect to the factor loading matrix `Λ`.
The gradient is output into `∇Q`.
The *∇Q* calculation is skipped if `∇Q = nothing`.

Returns the *Q* criterion value.
"""
criterion_and_gradient!

# fallback method that uses AutoDiff
function criterion_and_gradient!(∇Q, method::RotationMethod, Λ::AbstractMatrix)
    if !isnothing(∇Q)
        gradient!(Reverse, ∇Q, Base.Fix1(criterion, method), Λ)
    end
    return criterion(method, Λ)
end

"""
    isorthogonal(::RotationMethod)

Checks if the supplied rotation method is orthogonal.

## Examples
```jldoctest
julia> isorthogonal(Varimax())
true

julia> isorthogonal(Oblimax(orthogonal = false))
false
```
"""
isorthogonal(method::RotationMethod) = method isa RotationMethod{Orthogonal}

"""
    isoblique(::RotationMethod)

Checks if the supplied rotation method is oblique.

## Examples
```jldoctest
julia> isoblique(Varimax())
false

julia> isoblique(Oblimax(orthogonal = false))
true
```
"""
isoblique(method::RotationMethod) = method isa RotationMethod{Oblique}


"""
    rotation_type(::RotationMethod)

Return the rotation type for a given rotation method.

## Examples
```jldoctest
julia> rotation_type(Varimax())
Orthogonal

julia> rotation_type(Oblimin(gamma = 0.5))
Oblique
```
"""
rotation_type(::RotationMethod{RT}) where {RT} = RT

include("biquartimax.jl")
include("biquartimin.jl")
include("component_loss.jl")
include("crawford_ferguson.jl")
include("geomin.jl")
include("infomax.jl")
include("minimum_entropy.jl")
include("minimum_entropy_ratio.jl")
include("oblimax.jl")
include("target_rotation.jl")
include("oblimin.jl")
include("pattern_simplicity.jl")
include("quartimax.jl")
include("simplimax.jl")
include("tandem_criteria.jl")
include("varimax.jl")
