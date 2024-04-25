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

## Defining the rotation quality criterion
The easiest way to define a rotation method is to implement the
[`criterion_and_gradient!(::Nothing, ...)`](@ref criterion_and_gradient!) function,
which calculates the rotation quality criterion.
`nothing` fixed as the first argument specifies that
this implementation does not calculate the criterion gradient.

```jldoctest implementing_rotation_methods
julia> import FactorRotations: criterion_and_gradient!

julia> function criterion_and_gradient!(::Nothing, method::MyQuartimax, Λ::AbstractMatrix)
           return -sum(Λ .^ 4)
       end;

julia> criterion(MyQuartimax(), ones(10, 2))
-20.0
```

!!! note
    Since the algorithm in this package minimizes the criterion value, we have to make sure to return `-sum(...)` instead of the original criterion for Quartimax.

*FactorRotations.jl* will apply [Automatic Differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) to derive the gradient for your quality criterion
and use it during rotation optimization.

```jldoctest implementing_rotation_methods
julia> grad = fill(NaN, 10, 2);

julia> criterion_and_gradient!(grad, MyQuartimax(), ones(10, 2))
-20.0

julia> grad
10×2 Matrix{Float64}:
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0
 -4.0  -4.0

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
FactorRotation{Float64} with loading matrix:
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
FactorRotation{Float64} with loading matrix:
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
julia> isapprox(loadings(L_rotated), loadings(L_reference), atol = 1e-5)
true
```

## Defining the rotation quality gradient
When the gradient formula is available, [`criterion_and_gradient!`](@ref) method
could be modified to allow `∇Q::AbstractMatrix` as the first argument.
In this case *FactorRotations.jl* will expect that the `criterion_and_gradient!(∇Q, method, Λ)` call
sets `∇Q` to the ``∇Q(Λ)`` in-place and also returns ``Q(Λ)``.

This variant of `criterion_and_gradient!` allows reusing intermediate computations between
the criterion and its gradient, as well as providing more efficient gradient calculation
than the autodiff-based one.

Continuing the example of `MyQuartimax`:

```jldoctest implementing_rotation_methods
julia> import FactorRotations: criterion_and_gradient!

julia> function criterion_and_gradient!(∇Q::Union{AbstractMatrix, Nothing}, method::MyQuartimax, Λ::AbstractMatrix)
           Q = -sum(Λ .^ 4)
           if !isnothing(∇Q)
               ∇Q .= -Λ.^3
           end
           return Q
       end;
```

User-defined `criterion_and_gradient!` has priority over the default autodiff-based one:

```jldoctest implementing_rotation_methods
julia> grad = fill(NaN, 10, 2);

julia> criterion_and_gradient!(grad, MyQuartimax(), ones(10, 2))
-20.0

julia> grad
10×2 Matrix{Float64}:
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0
 -1.0  -1.0

```

[`rotate`](@ref) will now also use the custom [`criterion_and_gradient!`](@ref):

```jldoctest implementing_rotation_methods; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> L_rotated = rotate(L, MyQuartimax())
FactorRotation{Float64} with loading matrix:
8×2 Matrix{Float64}:
 0.898755  0.194824
 0.933943  0.129749
 0.902131  0.103865
 0.876508  0.171285
 0.315572  0.876476
 0.251123  0.773489
 0.198007  0.714678
 0.307857  0.659335

julia> isapprox(loadings(L_rotated), loadings(L_reference), atol = 1e-5)
true
```

Note that in our `criterion_and_gradient!` method gradient calculation is optional
and could be skipped by passing `nothing`. This variant of the method is used by [`criterion`](@ref):

```jldoctest implementing_rotation_methods; filter = r"([0-9]*)\.([0-9]{4})[0-9]+" => s"\1.\2"
julia> criterion(MyQuartimax(), L)
-2.8829327011730004
```
