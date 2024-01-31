"""
    FactorRotation{T <: Real}

A type holding results of a factor rotation.

## Fields
- `L`: The rotated factor loading matrix
- `T`: The factor rotation matrix
- `phi`: The factor correlation matrix
"""
struct FactorRotation{T}
    L::Matrix{T}
    T::Matrix{T}
    phi::Matrix{T}
end

function FactorRotation(L, T)
    return FactorRotation(L, T, T' * T)
end

function Base.show(io::IO, r::FactorRotation)
    println(io, "$(typeof(r)) with loading matrix:")
    show(io, "text/plain", r.L)
    return nothing
end

"""
    loadings(r::FactorRotation)

Return the rotated factor loading matrix from `r`.
"""
loadings(r::FactorRotation) = r.L

"""
    rotation(r::FactorRotation)

Return the factor rotation matrix from `r`.
"""
rotation(r::FactorRotation) = r.T

"""
    factor_correlation(r::FactorRotation)

Return the factor correlation matrix from `r`.
"""
factor_correlation(r::FactorRotation) = r.phi

"""
    rotate(Λ, method::RotationMethod; kwargs...)

Perform a rotation of the factor loading matrix `Λ` using a rotation `method`.

## Keyword arguments
- `atol`: Sets the absolute tolerance for convergence of the algorithm (default: 1e-6).
- `alpha`: Sets the inital value for alpha (default: 1).
- `maxiter1`: Controls the number of maximum iterations in the outer loop of the algorithm
              (default: 1000).
- `maxiter2`: Controls the number of maximum iterations in the inner loop of the algorithm
              (default: 10).
- `init`: A k-by-k matrix of starting values for the algorithm.
          If `init = nothing` (the default), the identity matrix will be used as starting
          values.
- `verbose`: Print logging statements (default: false)

## Examples
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> rotate(L, Varimax())
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
"""
function rotate(Λ, method; verbose = false, kwargs...)
    rotation = _rotate(Λ, method; verbose, kwargs...)

    if verbose
        @info "Rotation algorithm converged after $(length(rotation.iterations)) iterations." algorithm =
            typeof(method) criterion = last(rotation.iterations).Q
    end

    return FactorRotation(rotation.L, rotation.T)
end

function rotate(Λ, method::TandemCriteria; kwargs...)
    rotation_1 = rotate(Λ, TandemCriterionI(); kwargs...)
    reduced_loading_matrix = loadings(rotation_1)[:, 1:method.keep]
    rotation_2 = rotate(reduced_loading_matrix, TandemCriterionII(); kwargs...)
    return rotation_2
end

"""
    rotate!(Λ, method::RotationMethod; kwargs...)

Perform a rotation of the factor loading matrix `Λ` and overwrite `Λ` with the rotated
loading matrix.

For a list of available keyword arguments see [`rotate`](@ref).

## Examples
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2"
$(DEFINITION_L)
julia> rotate!(L, Quartimax())
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
"""
function rotate!(Λ, method; kwargs...)
    rot = rotate(Λ, method; kwargs...)
    Λ .= loadings(rot)
    return Λ
end

"""
    IterationState

A struct that holds the state of an iteration of the rotation algorithm.
"""
struct IterationState{T<:Real,V<:Real}
    alpha::T
    maxiter::Int
    Q::V
end

"""
    RotationState

A struct that holds the state of the rotation algorithm.

## Fields
- `init`: The initial rotation matrix
- `A`: The initial factor loading matrix
- `T`: The current state of the rotation matrix
- `L`: The current state of the rotated loading matrix
- `iterations`: A vector of [`IterationState`](@ref) that holds iteration states of the
                optimization.
"""
mutable struct RotationState{RT<:RotationType,T<:AbstractMatrix}
    init::T
    A::T
    T::T
    Ti::Union{Nothing,T}
    L::T
    iterations::Vector{IterationState}
end

function RotationState(T::Type{Orthogonal}, init, A)
    return RotationState{T,typeof(A)}(init, A, init, nothing, A * init, IterationState[])
end

function RotationState(T::Type{Oblique}, init, A)
    Ti = inv(init)
    return RotationState{T,typeof(A)}(init, A, init, Ti, A * Ti', IterationState[])
end

"""
    _rotate(A::AbstractMatrix, method::RotationMethod{Orthogonal}; kwargs...)

Implements the algorithm for factor rotation described in Bernaard & Jennrich (2005) for
orthogonal factor rotation.
"""
function _rotate(
    A::AbstractMatrix{TV},
    method::RotationMethod{RT};
    atol = 1e-6,
    alpha = 1,
    maxiter1 = 1000,
    maxiter2 = 10,
    init::Union{Nothing,AbstractMatrix} = nothing,
    verbose = false,
) where {RT,TV<:Real}
    state = initialize(RT, init, A)
    Q, ∇Q = criterion_and_gradient(method, state.L)

    # preallocate variables to avoid unnecessary allocations
    ft = Q
    G = gradient_f(state, ∇Q)
    Tt = similar(state.T)
    Gp = similar(state.T)
    s = zero(eltype(G))

    for i in 1:maxiter1
        project_G!(state, Gp, G)
        s = norm(Gp)

        isconverged(s, atol) && break

        alpha *= 2

        for _ in 1:maxiter2
            X = state.T - alpha * Gp
            Tt = project_X(state, X)
            update_state!(state, Tt)

            Q, ∇Q = criterion_and_gradient(method, state.L)

            if (Q < ft - 0.5 * s^2 * alpha)
                if verbose
                    @info "State at iteration $(i):" criterion = Q alpha = alpha
                end
                state.T = Tt
                ft = Q
                G = gradient_f(state, ∇Q)
                break
            else
                alpha /= 2
            end
        end

        iteration_state = IterationState(alpha, maxiter2, Q)
        push!(state.iterations, iteration_state)
    end

    isconverged(s, atol) || error("Algorithm did not converge after $(maxiter1) iterations")

    return state
end

"""
    initialize(init, A)

Initialize a [`RotationState`](@ref) with initial values `init` and original loading matrix
`A`. If `init = nothing`, the identity matrix will be used as initial values.
"""
function initialize(::Type{RT}, init, A::AbstractMatrix{TV}) where {RT<:RotationType,TV}
    _, k = size(A)

    if isnothing(init)
        T = Matrix{TV}(I, k, k)
    else
        T = init
    end

    if size(T) != (k, k)
        throw(ArgumentError("matrix of starting values must be of size ($k, $k)"))
    end

    return RotationState(RT, T, A)
end

function gradient_f(state::RotationState{Orthogonal}, ∇Q)
    @unpack A = state
    G = A' * ∇Q
    return G
end

function gradient_f(state::RotationState{Oblique}, ∇Q)
    @unpack L, Ti = state
    G = -(L' * ∇Q * Ti)'
    return G
end

"""
    project!(Gp, T, G)

Compute the projection `Gp` of `G` and store the results in `Gp`.
"""
function project_G!(state::RotationState{Orthogonal}, Gp, G)
    @unpack T = state
    M = T' * G
    S = (M + M') ./ 2
    Gp .= G - T * S
    return Gp
end

function project_G!(state::RotationState{Oblique}, Gp, G)
    @unpack T = state
    TG = T .* G
    Gp .= G - T * diagm(vec(sum(TG, dims = 1)))
    return Gp
end

"""
compute the projection `Tt` of `X`.
"""
function project_X(state::RotationState{Orthogonal}, X)
    @unpack U, Vt = svd(X)
    Tt = U * Vt
    return Tt
end

function project_X(state::RotationState{Oblique}, X)
    Xsq = X .^ 2
    v = 1 ./ sqrt.(sum(Xsq, dims = 1))
    Tt = X * diagm(vec(v))
    return Tt
end

"""
update the rotation state given a new projection `Tt`.
"""
function update_state!(state::RotationState{Orthogonal}, Tt)
    return mul!(state.L, state.A, Tt)
end

function update_state!(state::RotationState{Oblique}, Tt)
    state.Ti = inv(Tt)
    return mul!(state.L, state.A, state.Ti')
end

"""
    isconverged(s, atol)

determines the convergence status of the algorithm.
"""
isconverged(s, atol) = s < atol
