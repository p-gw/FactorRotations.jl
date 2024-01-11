# Implementing a rotation method

If you wish to implement your own factor rotation method or extend this package, you can do so in two ways:

- Implementing a rotation method without specifying the gradient
- Implementing a rotation method with gradients

In the following guide we will walk through both ways of implementing a rotation method.
As an example we will reimplement [`Quartimax`](@ref), which minimizes 

```math
Q = \sum_p \sum_k \lambda_{pk}^4
```

where ``p`` and ``k`` are the row and column indices of factor loading matrix ``\Lambda`` and ``\lambda_{pk}`` are the entries of the factor loading matrix. 

The first step to a custom implementation is to define a new `struct` for the rotation method.
FactorRotations.jl requires that all rotation methods inherit from [`RotationMethod`](@ref).
One must also specify whether the new method can be used for orthogonal rotation, oblique rotation, or both. 
For orthogonal rotation it is required that `T <: RotationMethod{Orthogonal}`.
Oblique rotations must satisfy `T <: RotationMethod{Oblique}`.
Methods that can be used for both orthogonal and oblique rotation are defined `T{RT} <: RotationMethod{RT}`.

Since Quartimax is an orthogonal rotation method, we define it as such.

```jldoctest implementing_rotation_methods
julia> using FactorRotations

julia> struct MyQuartimax <: RotationMethod{Orthogonal} end

```

## Gradient free methods
If no gradients are available, the easiest way to define a rotation method is to implement [`criterion`](@ref).
FactorRotations.jl will then use [Automatic Differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) to derive the gradients for your method.

```jldoctest implementing_rotation_methods
julia> import FactorRotations: criterion

julia> function criterion(method::MyQuartimax, Λ::AbstractMatrix)
           return -sum(Λ .^ 4)
       end;

julia> criterion(MyQuartimax(), ones(10, 2))
-20.0
```

!!! note
    Since the algorithm in this package minimizes the criterion value, we have to make sure to return `-sum(...)` instead of the original criterion for Quartimax. 

Once implemented, the required gradients for rotation will be automatically supplied by FactorRotations.jl.

```jldoctest implementing_rotation_methods
julia> criterion_and_gradient(MyQuartimax(), ones(10, 2))
(-20.0, [-4.0 -4.0; -4.0 -4.0; … ; -4.0 -4.0; -4.0 -4.0])
```

```jldoctest implementing_rotation_methods; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L = [
           0.830 -0.396
           0.818 -0.469
           0.777 -0.470
           0.798 -0.401
           0.786  0.500
           0.672  0.458
           0.594  0.444
           0.647  0.333
       ];

julia> L_rotated = rotate(L, MyQuartimax()) 
┌ Info: Rotation algorithm converged after 16 iterations.
│       algorithm: MyQuartimax
└       criterion: -4.090939184775285
8×2 Matrix{Float64}:
 0.898755  0.194824
 0.933943  0.129749
 0.902131  0.103864
 0.876508  0.171284
 0.315572  0.876476
 0.251123  0.773489
 0.198007  0.714678
 0.307857  0.659334
```

Checking against the [`Quartimax`](@ref) implementation shows that the results are approximately equal.

```jldoctest implementing_rotation_methods; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L_reference = rotate(L, Quartimax())
┌ Info: Rotation algorithm converged after 13 iterations.
│       algorithm: Quartimax
└       criterion: -1.0227347961934468
8×2 Matrix{Float64}:
 0.898755  0.194823
 0.933943  0.129748
 0.902132  0.103864
 0.876508  0.171284
 0.315572  0.876476
 0.251124  0.773489
 0.198008  0.714678
 0.307858  0.659334
```

```jldoctest implementing_rotation_methods
julia> isapprox(L_rotated, L_reference, atol = 1e-5)
true
```

## Methods with gradients
When gradients are available it can be helpful to implementing them directly.
In this case [`criterion_and_gradient`](@ref) can be directly implemented.
This can be beneficial for example if computation can be reused between the calculation of the criterion and its gradient.

Continuing the example of `MyQuartimax`:

```jldoctest implementing_rotation_methods
julia> import FactorRotations: criterion_and_gradient

julia> function criterion_and_gradient(method::MyQuartimax, Λ::AbstractMatrix)
           Q = criterion(method, Λ)
           ∇Q = -Λ.^3
           return Q, ∇Q
       end;
```

Evaluating the function will now use the newly defined method, 

```jldoctest implementing_rotation_methods
julia> criterion_and_gradient(MyQuartimax(), ones(10, 2))
(-20.0, [-1.0 -1.0; -1.0 -1.0; … ; -1.0 -1.0; -1.0 -1.0])
```

Again, this method can be simply used with [`rotate`](@ref), now using the custom [`criterion_and_gradient`](@ref).

```jldoctest implementing_rotation_methods; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L_rotated = rotate(L, MyQuartimax())
┌ Info: Rotation algorithm converged after 12 iterations.
│       algorithm: MyQuartimax
└       criterion: -4.090939184774669
8×2 Matrix{Float64}:
 0.898755  0.194824
 0.933943  0.129749
 0.902131  0.103865
 0.876508  0.171285
 0.315572  0.876476
 0.251123  0.773489
 0.198007  0.714678
 0.307857  0.659335

julia> isapprox(L_rotated, L_reference, atol = 1e-5)
true
```
