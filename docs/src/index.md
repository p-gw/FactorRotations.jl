# FactorRotations.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://p-gw.github.io/FactorRotations.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://p-gw.github.io/FactorRotations.jl/dev/)
[![Build Status](https://github.com/p-gw/FactorRotations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/p-gw/FactorRotations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/p-gw/FactorRotations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/p-gw/FactorRotations.jl)

[FactorRotations.jl](https://github.com/p-gw/FactorRotations.jl) implements factor rotations by the gradient projections algorithms described
in Bernaards & Jennrich (2005).

## Installation
To install FactorRotations.jl you can use the Julia package manager,

```julia
] add FactorRotations
```

## Getting started
*FactorRotations.jl* provides methods to rotate factor loading matrices, e.g. from
exploratory factor analysis or principle component analysis.

Assume you aquired a factor loading matrix `L` then you can rotate the matrix by calling
the `rotate` function. The `rotate` function takes the factor loading matrix as the first
argument and an instance of a rotation method as the second argument.

```julia
L = [
    0.830 -0.396
    0.818 -0.469
    0.777 -0.470
    0.798 -0.401
    0.786  0.500
    0.672  0.458
    0.594  0.444
    0.647  0.333
]

rotate(L, Varimax())
```

For a complete list of available methods see the [Rotation Methods](@ref rotation_methods) section of the documentation.

For a fully worked example see the [Guides](@ref guides) section of the documentation.

# References
Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms and software for arbitrary rotation criteria in factor analysis. *Educational and psychological measurement, 65*(5), 676-696.
