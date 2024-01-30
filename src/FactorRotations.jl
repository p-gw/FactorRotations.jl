module FactorRotations

using Enzyme
using FillArrays
using LinearAlgebra
using LogExpFunctions
using SimpleUnPack
using Statistics

import LinearAlgebra: rotate!

export FactorRotation, loadings, rotation, factor_correlation
export rotate, rotate!
export criterion, criterion_and_gradient

export kaiser_normalize, kaiser_denormalize
export kaiser_normalize!, kaiser_denormalize!

export RotationMethod
export Orthogonal, Oblique
export isorthogonal, isoblique

export Biquartimax
export ComponentLoss, KatzRohlf, LinearRightConstant, Concave, Absolmin
export CrawfordFerguson
export Cubimax
export Geomin
export PatternSimplicity
export Infomax
export MinimumEntropy
export MinimumEntropyRatio
export Oblimax
export Oblimin
export Quartimax
export Simplimax
export TandemCriterionI, TandemCriterionII, TandemCriteria
export TargetRotation
export Varimax

include("utils.jl")
include("normalize.jl")
include("rotation_types.jl")
include("methods/methods.jl")
include("rotate.jl")

end
