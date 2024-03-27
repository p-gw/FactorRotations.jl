```@meta
CurrentModule = FactorRotations
```

# API

## Rotation criteria

```@docs
Absolmin
Biquartimax
ComponentLoss
Concave
CrawfordFerguson
Infomax
KatzRohlf
LinearRightConstant
MinimumEntropy
MinimumEntropyRatio
Oblimax
Oblimin
PatternSimplicity
Quartimax
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
```
