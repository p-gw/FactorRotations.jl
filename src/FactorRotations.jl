module FactorRotations

using Folds
using FillArrays
using LinearAlgebra
using LogExpFunctions
using Logging
using SimpleUnPack
using Statistics

import LinearAlgebra: rotate!, reflect!

export setverbosity!

export FactorRotation, loadings, rotation, factor_correlation
export rotate, rotate!
export criterion, criterion_and_gradient!

export reflect, reflect!

export kaiser_normalize, kaiser_denormalize
export kaiser_normalize!, kaiser_denormalize!

export RotationMethod
export Orthogonal, Oblique
export isorthogonal, isoblique
export rotation_type

export Biquartimax
export Biquartimin
export ComponentLoss, KatzRohlf, LinearRightConstant, Concave, Absolmin
export CrawfordFerguson
export Cubimax
export Equamax
export Geomin
export Parsimax
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

const VERBOSITY = Ref(false)

"""
    setverbosity!(::Bool)

Sets the global verbosity level of the package.
If set to `false` (the default), package functions will not log `@info` statements.
If set to `true`, package functions will provide `@info` statements.
"""
function setverbosity!(verbose::Bool)
    @info "Logging is $(verbose ? "enabled" : "disabled") globally."
    VERBOSITY[] = verbose
    return nothing
end

struct AutodiffBackend{B}
    AutodiffBackend(backend::Symbol) = new{backend}()
end

const AUTODIFF_BACKEND = Ref{AutodiffBackend}(AutodiffBackend(:Enzyme))

"""
    set_autodiff_backend(backend::Symbol)

Sets the *automatic differentiation* backend.

Automatic differentiation is used by the fallback `criterion_and_gradient!()` implementation.
Currently, only `:Enzyme` backend is supported.

Note that to actually enable the differentiation,
the corresponding autodiff package must be loaded first
(e.g. `using Enzyme`)
"""
function set_autodiff_backend(backend::Symbol)
    @info "$(@__MODULE__) autodiff backend set to $(backend)."
    AUTODIFF_BACKEND[] = AutodiffBackend(backend)
    return nothing
end

include("utils.jl")
include("normalize.jl")
include("rotation_types.jl")
include("methods/methods.jl")
include("rotate.jl")
include("reflect.jl")

end
