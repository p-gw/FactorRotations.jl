# Rotate an existing loading matrix

This guide considers the basic use case of FactorRotations.jl: Given an existing factor loading matrix `L`, calculate the rotation of the loading matrix according to some rotation criterion. In this example we will first consider the simple case of orthogonal _Varimax_ rotation. At a later stage we will see how to easily switch the factor rotation criterion to arrive at a different rotation.

First, we assume a factor loading matrix `L`. 
In this example we will use the loading matrix given by Bernaard & Jennrich (2005),

```jldoctest basic_example
julia> using FactorRotations

julia> L = [
           0.830 -0.396
           0.818 -0.469
           0.777 -0.470
           0.798 -0.401
           0.786  0.500
           0.672  0.458
           0.594  0.444
           0.647  0.333
       ]
8×2 Matrix{Float64}:
 0.83   -0.396
 0.818  -0.469
 0.777  -0.47
 0.798  -0.401
 0.786   0.5
 0.672   0.458
 0.594   0.444
 0.647   0.333
```

Rotating the loading matrix consists of a single call to [`rotate`](@ref). This function takes the unrotated loading matrix as the first argument, and an instance of a factor rotation method as a second argument.

For clarity we first set up our [`Varimax`](@ref) rotation method,

```jldoctest basic_example
julia> criterion = Varimax()
Varimax()
```

Finally we perform the rotation using [`rotate`](@ref),

```jldoctest basic_example; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L_rotated = rotate(L, criterion)
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.886061  0.246196
 0.924934  0.183253
 0.894664  0.155581
 0.865205  0.221416
 0.264636  0.893176
 0.206218  0.786653
 0.156572  0.724884
 0.269424  0.67595
```

Different rotation can be achieved by simply changing `criterion` or passing it directly to [`rotate`](@ref).

```jldoctest basic_example; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L_rotated = rotate(L, MinimumEntropy())
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.90711   0.151221
 0.939117  0.084524
 0.906093  0.0602051
 0.883753  0.128783
 0.357504  0.860225
 0.28816   0.760468
 0.232268  0.704289
 0.339319  0.643709 
```

## In-place rotation
In some cases it can be useful to modify `L` directly. 
For this use case the package provides an in-place rotation, [`rotate!`](@ref) with the same function signature as before.

```jldoctest basic_example; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> rotate!(L, MinimumEntropy())
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.90711   0.151221
 0.939117  0.084524
 0.906093  0.0602051
 0.883753  0.128783
 0.357504  0.860225
 0.28816   0.760468
 0.232268  0.704289
 0.339319  0.643709 

julia> L == L_rotated
true
```
