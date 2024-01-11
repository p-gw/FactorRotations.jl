const DEFINITION_L = """
julia> L = [
           0.830 -0.396
           0.818 -0.469
           0.777 -0.470
           0.798 -0.401
           0.786  0.500
           0.672  0.458
           0.594  0.444
           0.647  0.333
       ];
"""

"""
    RotationMethod{T<:RotationType}

An abstract type representing a factor rotation method.

Each implementation of `M <: RotationMethod` must implement at least one of the following methods:
- [`criterion`](@ref)
- [`criterion_and_gradient`](@ref)

If [`criterion`](@ref) is implemented, gradients are calculated by automatic differentiation.
"""
abstract type RotationMethod{T<:RotationType} end

"""
    criterion(method::RotationMethod, Λ::Abstractmatrix{<:Real})

Calculate the criterion of a given `method` with respect to the factor loading matrix `Λ`.
"""
function criterion end

"""
    criterion_and_gradient(method::RotationMethod, Λ::AbstractMatrix{<:Real})

Calculate the criterion and gradient of a given `method` with respect to the factor loading
matrix `Λ`.

Returns a Tuple with the criterion value as the first element and gradient as the second
element.
"""
function criterion_and_gradient(method::RotationMethod, Λ::AbstractMatrix)
    Q = criterion(method, Λ)
    ∇Q = gradient(Reverse, x -> criterion(method, x), Λ)
    return Q, ∇Q
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

include("biquartimax.jl")
include("crawford_ferguson.jl")
include("cubimax.jl")
include("infomax.jl")
include("minimum_entropy.jl")
include("minimum_entropy_ratio.jl")
include("oblimax.jl")
include("target_rotation.jl")
include("oblimin.jl")
include("quartimax.jl")
include("simplimax.jl")
include("varimax.jl")
