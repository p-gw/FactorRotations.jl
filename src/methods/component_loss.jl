abstract type AbstractComponentLoss{RT} <: RotationMethod{RT} end

criterion(method::AbstractComponentLoss, Λ::AbstractMatrix) = -sum(method.loss, Λ)

struct ComponentLoss{RT,F} <: AbstractComponentLoss{RT}
    loss::F
    function ComponentLoss(loss::F; orthogonal = false) where {F}
        T = orthogonal ? Orthogonal : Oblique
        return new{T,F}(loss)
    end
end

struct KatzRohlf{F} <: AbstractComponentLoss{Orthogonal}
    bandwidth::Float64
    loss::F
    function KatzRohlf(bandwidth)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = 1 - exp(-(x / bandwidth)^2)
        return new{typeof(loss)}(bandwidth, loss)
    end
end

struct LinearRightConstant{F} <: AbstractComponentLoss{Orthogonal}
    bandwidth::Float64
    loss::F
    function LinearRightConstant(bandwidth)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = abs(x) > bandwidth ? 1.0 : (x / bandwidth)^2
        return new{typeof(loss)}(bandwidth, loss)
    end
end

struct Concave{F} <: AbstractComponentLoss{Oblique}
    bandwidth::Float64
    loss::F
    function Concave(bandwidth = 1)
        bandwidth > 0 || throw(ArgumentError("bandwidth must be positive"))
        loss(x) = 1 - exp(-abs(x) / bandwidth)
        return new{typeof(loss)}(bandwidth, loss)
    end
end

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
