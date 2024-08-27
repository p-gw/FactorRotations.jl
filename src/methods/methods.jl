"""
    RotationMethod{T<:RotationType}

An abstract type representing a factor rotation method.

Each implementation of `M <: RotationMethod` must provide [`criterion_and_gradient!`](@ref)
method.
"""
abstract type RotationMethod{T<:RotationType} end

"""
    criterion(method::RotationMethod, Λ::Abstractmatrix{<:Real})

Calculate the criterion of a given `method` with respect to the factor loading matrix `Λ`.

The method is just a wrapper for a [`criterion_and_gradient!(nothing, method, Λ)`](@ref criterion_and_gradient!) call.
"""
criterion(method::RotationMethod, Λ::AbstractMatrix{<:Real}) = criterion_and_gradient!(nothing, method, Λ)

"""
    criterion_and_gradient!(∇Q::Union{AbstractMatrix{<:Real}, Nothing},
                            method::RotationMethod, Λ::AbstractMatrix{<:Real})

Calculate the quality criterion *Q* and its gradient for a given `method`
with respect to the factor loading matrix `Λ`.
The gradient is output into `∇Q` matrix, which should have the same dimensions as `Λ`.
The `∇Q` calculation is skipped if `∇Q ≡ nothing`.

Returns the *Q* criterion value.
"""
criterion_and_gradient!

OptionalGradient = Union{AbstractMatrix, Nothing}

# fallback method when autodiff backend is not available
function autodiff_gradient!(_::AutodiffBackend{B}, ∇Q::AbstractMatrix, method::RotationMethod, Λ::AbstractMatrix) where B
    if B == :Enzyme
        error("Enzyme.jl autodiff backend is not loaded. Have you run \"using Enzyme\"?")
    else
        error("$(B) autodiff backend is not supported.")
    end
end

autodiff_gradient!(∇Q::AbstractMatrix, method::RotationMethod, Λ::AbstractMatrix) =
    autodiff_gradient!(AUTODIFF_BACKEND[], ∇Q, method, Λ)

# fallback method that applies auto-diff to criterion() call
function criterion_and_gradient!(∇Q::OptionalGradient, method::RotationMethod, Λ::AbstractMatrix)
    if !isnothing(∇Q)
        autodiff_gradient!(∇Q, method, Λ)
    else
        error("$(typeof(method)) does not implement neither criterion_and_gradient!(∇Q, ...) nor criterion_and_gradient!(nothing, ...) methods.")
    end
    return criterion(method, Λ)
end

"""
    isorthogonal(::RotationMethod)

Checks if the supplied rotation method is orthogonal.

## Examples
```jldoctest
julia> isorthogonal(Varimax())
true

julia> isorthogonal(Oblimax(orthogonal = false))
false
```
"""
isorthogonal(method::RotationMethod) = method isa RotationMethod{Orthogonal}

"""
    isoblique(::RotationMethod)

Checks if the supplied rotation method is oblique.

## Examples
```jldoctest
julia> isoblique(Varimax())
false

julia> isoblique(Oblimax(orthogonal = false))
true
```
"""
isoblique(method::RotationMethod) = method isa RotationMethod{Oblique}


"""
    rotation_type(::RotationMethod)

Return the rotation type for a given rotation method.

## Examples
```jldoctest
julia> rotation_type(Varimax())
Orthogonal

julia> rotation_type(Oblimin(gamma = 0.5))
Oblique
```
"""
rotation_type(::RotationMethod{RT}) where {RT} = RT

"""
    weighted_sums_criterion_and_gradient!(
        [∇Q::AbstractMatrix],
        Λ::AbstractMatrix{<:Real},
        columns_weight::Number, rows_weight::Number
    )

Calculate the value and the gradient of the criterion, which
is based on the weighted column- and row-wise sums of *Λ²*.

Specifically,
```math
Q(Λ) = \\frac{c}{4} ∑ⱼᵐ \\left(∑ᵢⁿ Λ²_{i,j}\\right)² +
       \\frac{r}{4} ∑ᵢⁿ \\left(∑ⱼᵐ Λ²_{i,j}\\right)² -
       \\frac{1}{4} ∑ᵢⁿ∑ⱼᵐ Λ⁴,
```
where *c* is `columns_weight` and *r* is `rows_weight`.

The gradient is output into `∇Q` matrix, which should have the same dimensions as `Λ`.
The `∇Q` calculation is skipped if `∇Q ≡ nothing`.

This function is used by multiple rotation methods, such as [`CrawfordFerguson`](@ref),
[`Equamax`](@ref), [`Oblimin`](@ref), and [`Parsimax`](@ref).
"""
function weighted_sums_criterion_and_gradient!(
    ∇Q::Union{Nothing, AbstractMatrix},
    Λ::AbstractMatrix{<:Real},
    columns_weight::Number, rows_weight::Number
)
    Λsq = !isnothing(∇Q) ? ∇Q : similar(Λ)
    Λsq .= Λ .^ 2

    Λsq_rowsum = sum(Λsq, dims=1)
    Λsq_colsum = sum(Λsq, dims=2)

    Q = (columns_weight * sum(abs2, Λsq_colsum) + rows_weight * sum(abs2, Λsq_rowsum) - sum(abs2, Λsq)) / 4
    if !isnothing(∇Q)
        # ∇Q === Λsq
        # weighted Λ² columns and rows sum at each position - Λ²
        ∇Q .= (columns_weight .* Λsq_colsum) .+
              (rows_weight .* Λsq_rowsum) .- Λsq
        ∇Q .*= Λ
    end
    return Q
end

include("biquartimax.jl")
include("biquartimin.jl")
include("component_loss.jl")
include("crawford_ferguson.jl")
include("equamax.jl")
include("geomin.jl")
include("infomax.jl")
include("minimum_entropy.jl")
include("minimum_entropy_ratio.jl")
include("oblimax.jl")
include("target_rotation.jl")
include("oblimin.jl")
include("parsimax.jl")
include("pattern_simplicity.jl")
include("quartimax.jl")
include("simplimax.jl")
include("tandem_criteria.jl")
include("varimax.jl")
