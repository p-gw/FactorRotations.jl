# Rotation Criteria

FactorRotations.jl implements both orthogonal and oblique rotation criteria. The following tables list the available methods and their respective references.

## Orthogonal rotation criteria

criterium                     | reference                  | note
----------------------------- | -------------------------- | -------------------------------------------------------
[`Biquartimax`](@ref)         |                            | equivalent to `Oblimin(gamma = 0.5, orthogonal = true)`
[`CrawfordFerguson`](@ref)    | Crawford & Ferguson (1970) |
[`Infomax`](@ref)             | Browne (2001)              | based on the unpublished manuscript McKeon (1968)
[`MinimumEntropyRatio`](@ref) |                            |
[`MinimumEntropy`](@ref)      |                            |
[`Oblimax`](@ref)             |                            |
[`Oblimin`](@ref)             |                            | equivalent to `Oblimin` for orthogonal rotation
[`Quartimax`](@ref)           | Neuhaus & Wrigley (1954)   | equivalent to `Oblimin(gamma = 0, orthogonal = true)`
[`TargetRotation`](@ref)      |                            |
[`Varimax`](@ref)             | Kaiser (1958)              | equivalent to `Oblimin(gamma = 1, orthogonal = true)`

## References

Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. _Multivariate behavioral research, 36_(1), 111-150.

Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal rotation. _Psychometrika, 35_, 321-332.

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. _Psychometrika, 23_(3), 187-200.

Neuhaus, J. O., & Wrigley, C. (1954). The quartimax method: An analytic approach to orthogonal simple structure. _British Journal of Statistical Psychology, 7_(2), 81-91.
