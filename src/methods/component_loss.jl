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
        loss(x) = 1 - exp(-(x / bandwidth)^2)
        return new{typeof(loss)}(bandwidth, loss)
    end
end

struct LinearRightConstant{F} <: AbstractComponentLoss{Orthogonal}
    bandwidth::Float64
    loss::F
    function LinearRightConstant(bandwidth)
        loss(x) = abs(x) <= bandwidth ? (x / bandwidth)^2 : 1
        return new{typeof(loss)}(bandwidth, loss)
    end
end
