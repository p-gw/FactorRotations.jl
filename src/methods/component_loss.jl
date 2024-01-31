"""
    AbstractComponentLoss{RT} <: RotationMethod{RT}

An abstract type representing component loss functions.

Implementing a custom component loss `T` requires that `T <: AbstractComponentLoss`.
"""
abstract type AbstractComponentLoss{RT} <: RotationMethod{RT} end

criterion(method::AbstractComponentLoss, Λ::AbstractMatrix) = -sum(method.loss, Λ)

"""
    ComponentLoss(loss::Function; orthogonal = false)

A generic implementation of the component loss factor rotation method.
`loss` defines the loss function that is applied to the components of the loading matrix.

## Keyword arguments
- `orthogonal`: If `orthogonal = true` an orthogonal rotation is performed, an oblique
   rotation otherwise. (default: `false`)

## Details
The component loss factor rotation applies a loss function to each element of the factor
loading matrix. Then the following criterion is minimized:

```math
Q(\\Lambda) = \\sum_i \\sum_j h(\\lambda_{ij})
```

## Examples
### Quartimax as a component loss
```jldoctest; filter = [r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2", r"var\\"(.*)\\"" => s""]
$(DEFINITION_L)
julia> quartimax_loss = ComponentLoss(x -> x^4, orthogonal = true);

julia> L_component_loss = rotate(L, quartimax_loss)
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

julia> L_quartimax = rotate(L, Quartimax())
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

julia> isapprox(loadings(L_component_loss), loadings(L_quartimax), atol = 1e-5)
true
```
"""
struct ComponentLoss{RT,F} <: AbstractComponentLoss{RT}
    loss::F
    function ComponentLoss(loss::F; orthogonal = false) where {F}
        T = orthogonal ? Orthogonal : Oblique
        return new{T,F}(loss)
    end
end

"""
    KatzRohlf(bandwidth)

A component loss criterion with loss function

```math
h(\\lambda) = 1 - \\exp(-(\\frac{\\lambda}{b})^2)
```

where ``b`` is the bandwidth parameter.
"""
struct KatzRohlf{F} <: AbstractComponentLoss{Orthogonal}
    bandwidth::Float64
    loss::F
    function KatzRohlf(bandwidth)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = 1 - exp(-(x / bandwidth)^2)
        return new{typeof(loss)}(bandwidth, loss)
    end
end

"""
    LinearRightConstant(bandwidth)

The linear right constant component loss factor rotation criterion.
It has the loss function

```math
h(\\lambda) = \\begin{cases}
    (\\frac{\\lambda}{b})^2&\\text{if } |\\lambda| \\leq b \\\\
    1 &\\text{if } \\lambda > b
\\end{cases}
```

where ``b`` is the bandwidth parameter.
"""
struct LinearRightConstant{F} <: AbstractComponentLoss{Orthogonal}
    bandwidth::Float64
    loss::F
    function LinearRightConstant(bandwidth)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = abs(x) > bandwidth ? 1.0 : (x / bandwidth)^2
        return new{typeof(loss)}(bandwidth, loss)
    end
end

"""
    Concave(bandwidth = 1)

The simple concave component loss factor rotation criterion.
It has the loss function

```math
h(\\lambda) = 1 - \\exp(-\\frac{|\\lambda|}{b})
```

where ``b`` is the bandwidth parameter.
"""
struct Concave{F} <: AbstractComponentLoss{Oblique}
    bandwidth::Float64
    loss::F
    function Concave(bandwidth = 1)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = 1 - exp(-abs(x) / bandwidth)
        return new{typeof(loss)}(bandwidth, loss)
    end
end

"""
    Absolmin(epsilon = 0)

The Absolmin component loss factor rotation criterion.
It has the loss function

```math
h(\\lambda) = |\\lambda|
```
"""
struct Absolmin{F} <: AbstractComponentLoss{Oblique}
    epsilon::Float64
    loss::F
    function Absolmin(epsilon = 0)
        epsilon >= 0 || throw(ArgumentError("epsilon must be non-negative"))
        b = 1 / (2 * epsilon)
        a = epsilon - b * epsilon^2
        loss(x) = abs(x) > epsilon ? abs(x) : a + b * abs2(x)
        return new{typeof(loss)}(epsilon, loss)
    end
end
