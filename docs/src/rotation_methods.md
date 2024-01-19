# Rotation Criteria

FactorRotations.jl implements both orthogonal and oblique rotation criteria. The following tables list the available methods and their respective references.

## Orthogonal rotation criteria

criterium                     | reference              | note
----------------------------- | ---------------------- | -------------------------------------------------------
[`Biquartimax`](@ref)         |                        | equivalent to `Oblimin(gamma = 0.5, orthogonal = true)`
[`CrawfordFerguson`](@ref)    | [crawford1970](@citet) |
[`Infomax`](@ref)             | [browne2001](@citet)   | based on the unpublished manuscript McKeon (1968)
[`KatzRohlf`](@ref)           |                        |
[`LinearRightConstant`](@ref) | [jennrich2004](@citet) |
[`MinimumEntropyRatio`](@ref) | [mccammon1966](@citet) |
[`MinimumEntropy`](@ref)      | [jennrich2004](@citet) |
[`Oblimax`](@ref)             |                        |
[`Oblimin`](@ref)             |                        | equivalent to `Oblimin` for orthogonal rotation
[`Quartimax`](@ref)           | [neuhaus1954](@citet)  | equivalent to `Oblimin(gamma = 0, orthogonal = true)`
[`TandemCriteria`](@ref)      | [comrey1967](@citet)   |
[`TandemCriterionII`](@ref)   | [comrey1967](@citet)   |
[`TandemCriterionI`](@ref)    | [comrey1967](@citet)   |
[`TargetRotation`](@ref)      |                        |
[`Varimax`](@ref)             | [kaiser1958](@citet)   | equivalent to `Oblimin(gamma = 1, orthogonal = true)`

## Oblique rotation criteria

criterium                   | reference                            | note
--------------------------- | ------------------------------------ | ---------------------------
[`Absolmin`](@ref)          | [jennrich2006](@citet)               |
[`ComponentLoss`](@ref)     | [jennrich2004, jennrich2006](@citet) | both orthogonal and oblique
[`Concave`](@ref)           | [jennrich2006](@citet)               |
[`PatternSimplicity`](@ref) | [bentler1977](@citet)                | both orthogonal and oblique

## References

```@bibliography
Pages = ["rotation_methods.md"]
```
