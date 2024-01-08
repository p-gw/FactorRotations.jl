module FactorRotations

using FillArrays
using LinearAlgebra
using LogExpFunctions
using SimpleUnPack
using Statistics

import LinearAlgebra: rotate!

export rotate, rotate!
export criterion_and_gradient

export Orthogonal, Oblique
export isorthogonal, isoblique

export Biquartimax
export CrawfordFerguson
export Infomax
export MinimumEntropy
export MinimumEntropyRatio
export Oblimax
export Oblimin
export Quartimax
export TargetRotation
export Varimax

include("utils.jl")
include("rotation_types.jl")
include("methods/methods.jl")
include("rotate.jl")

end
