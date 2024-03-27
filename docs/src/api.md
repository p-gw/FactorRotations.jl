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
Geomin
Infomax
KatzRohlf
LinearRightConstant
MinimumEntropy
MinimumEntropyRatio
Oblimax
Oblimin
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
rotate
rotate!
rotation
rotation_type
setverbosity!
```

## Internals

```@docs
FactorRotation
Oblique
Orthogonal
RotationMethod
criterion
criterion_and_gradient
ConvergenceError
```
