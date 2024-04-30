# [Rotation Methods](@id rotation_methods)

FactorRotations.jl implements both orthogonal and oblique rotation criteria. The following tables list the available methods and their respective references.

## Orthogonal rotation criteria

| criterium                     | reference                            | note                                                    |
| ----------------------------- | ------------------------------------ | ------------------------------------------------------- |
| [`Biquartimax`](@ref)         |                                      | equivalent to `Oblimin(gamma = 0.5, orthogonal = true)` |
| [`Biquartimin`](@ref)         | [jennrich2011](@citet)               |
| [`ComponentLoss`](@ref)       | [jennrich2004, jennrich2006](@citet) |
| [`CrawfordFerguson`](@ref)    | [crawford1970](@citet)               |
| [`Infomax`](@ref)             | [browne2001](@citet)                 | based on the unpublished manuscript McKeon (1968)       |
| [`KatzRohlf`](@ref)           |                                      |
| [`LinearRightConstant`](@ref) | [jennrich2004](@citet)               |
| [`MinimumEntropyRatio`](@ref) | [mccammon1966](@citet)               |
| [`MinimumEntropy`](@ref)      | [jennrich2004](@citet)               |
| [`Oblimax`](@ref)             |                                      |
| [`Oblimin`](@ref)             |                                      | equivalent to `Oblimin` for orthogonal rotation         |
| [`PatternSimplicity`](@ref)   | [bentler1977](@citet)                |
| [`Quartimax`](@ref)           | [neuhaus1954](@citet)                | equivalent to `Oblimin(gamma = 0, orthogonal = true)`   |
| [`TandemCriteria`](@ref)      | [comrey1967](@citet)                 |
| [`TandemCriterionII`](@ref)   | [comrey1967](@citet)                 | second step of [`TandemCriteria`](@ref)                 |
| [`TandemCriterionI`](@ref)    | [comrey1967](@citet)                 | first step of [`TandemCriteria`](@ref)                  |
| [`TargetRotation`](@ref)      |                                      |
| [`Varimax`](@ref)             | [kaiser1958](@citet)                 | equivalent to `Oblimin(gamma = 1, orthogonal = true)`   |

## Oblique rotation criteria

| criterium                   | reference                            | note                                              |
| --------------------------- | ------------------------------------ | ------------------------------------------------- |
| [`Absolmin`](@ref)          | [jennrich2006](@citet)               |
| [`Biquartimin`](@ref)       | [jennrich2011](@citet)               |
| [`ComponentLoss`](@ref)     | [jennrich2004, jennrich2006](@citet) |
| [`Concave`](@ref)           | [jennrich2006](@citet)               |
| [`CrawfordFerguson`](@ref)  | [crawford1970](@citet)               |
| [`Geomin`](@ref)            |                                      |
| [`Infomax`](@ref)           | [browne2001](@citet)                 | based on the unpublished manuscript McKeon (1968) |
| [`Oblimax`](@ref)           |                                      |
| [`Oblimin`](@ref)           |                                      |
| [`PatternSimplicity`](@ref) | [bentler1977](@citet)                |
| [`Simplimax`](@ref)         |                                      |
| [`TargetRotation`](@ref)    |                                      |

## References

```@bibliography
Pages = ["rotation_methods.md"]
```
