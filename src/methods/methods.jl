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

Each implementation of M <: RotationMethod must implement [`criterion_and_gradient`](@ref).
"""
abstract type RotationMethod{T<:RotationType} end

"""
    criterion_and_gradient(method::RotationMethod, Λ::AbstractMatrix{<:Real})

Calculate the criterion and gradient of a given `method` with respect to the factor loading
matrix `Λ`.

Returns a Tuple with the criterion value as the first element and gradient as the second
element.
"""
function criterion_and_gradient end

isorthogonal(method::RotationMethod) = method isa RotationMethod{Orthogonal}
isoblique(method::RotationMethod) = method isa RotationMethod{Oblique}

include("biquartimax.jl")
include("crawford_ferguson.jl")
include("infomax.jl")
include("minimum_entropy.jl")
include("minimum_entropy_ratio.jl")
include("oblimax.jl")
include("target_rotation.jl")
include("oblimin.jl")
include("quartimax.jl")
include("varimax.jl")
