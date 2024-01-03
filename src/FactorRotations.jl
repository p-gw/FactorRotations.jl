module FactorRotations

using FillArrays
using LinearAlgebra
using LogExpFunctions
using SimpleUnPack

export rotate, rotate!
export criterion_and_gradient

export Orthogonal, Oblique
export isorthogonal, isoblique

export Biquartimax
export CrawfordFerguson
export Oblimin
export Quartimax
export Varimax
export Infomax

include("utils.jl")
include("rotation_types.jl")
include("methods/methods.jl")
include("rotate.jl")

end
