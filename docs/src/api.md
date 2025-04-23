```@meta
CurrentModule = FactorRotations
```

# API

## Rotation criteria

Each of these rotation methods defines its specific quality criterion
``Q(Λ)``, where ``Λ = L⋅R ∈ ℝ^{p \times k}`` is the rotated factor loading matrix,
``L ∈ ℝ^{p \times k}`` is the original loading matrix, and ``R ∈ ℝ^{k \times k}``
is the *rotation matrix*.
The goal is to find the optimal rotation matrix ``R`` (either [*orthogonal*](@ref rotation_orthogonal)
or [*oblique*](@ref rotation_oblique)) that will minimize ``Q(Λ)``.

### Component loss-based

```@docs
AbstractComponentLoss
ComponentLoss
Absolmin
Concave
KatzRohlf
LinearRightConstant
```

### Based on column and row variances

These rotation methods optimize the variance of columns and/or rows of
the squared loadings matrix. The methods differ by how much weight is
assigned to the column and row variances.

```@docs
CrawfordFerguson
Oblimin
Biquartimax
Biquartimin
Equamax
Oblimax
Parsimax
Quartimax
Varimax
```

### Other

```@docs
Geomin
Infomax
MinimumEntropy
MinimumEntropyRatio
PatternSimplicity
Simplimax
TandemCriteria
TandemCriterionI
TandemCriterionII
TargetRotation
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
