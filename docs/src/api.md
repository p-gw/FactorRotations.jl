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
setverbosity!
rotate
rotate!
isorthogonal
isoblique
rotation_type
loadings
rotation
factor_correlation
```

## Internals

```@docs
FactorRotation
RotationMethod
Orthogonal
Oblique
criterion
criterion_and_gradient
ConvergenceError
```
