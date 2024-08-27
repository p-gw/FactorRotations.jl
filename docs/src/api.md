```@meta
CurrentModule = FactorRotations
```

# API

## Rotation criteria

```@docs
Absolmin
Biquartimax
Biquartimin
ComponentLoss
Concave
CrawfordFerguson
Equamax
Geomin
Infomax
KatzRohlf
LinearRightConstant
MinimumEntropy
MinimumEntropyRatio
Oblimax
Oblimin
Parsimax
PatternSimplicity
Quartimax
Simplimax
TandemCriteria
TandemCriterionI
TandemCriterionII
TargetRotation
Varimax
```

## User Interface

```@docs
factor_correlation
isoblique
isorthogonal
kaiser_denormalize
kaiser_denormalize!
kaiser_normalize
kaiser_normalize!
loadings
reflect
reflect!
rotate
rotate!
rotation
rotation_type
setverbosity!
FactorRotations.set_autodiff_backend
```

## Internals

```@docs
FactorRotation
Oblique
Orthogonal
RotationMethod
criterion
criterion_and_gradient!
```
