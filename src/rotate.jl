"""
    rotate(Λ, method::RotationMethod; kwargs...)
    rotate!(Λ, method::RotationMethod; kwargs...)

Perform a rotation of the factor loading matrix `Λ`.

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
"""
function rotate(Λ, method; kwargs...)
    rotation = _rotate(Λ, method; kwargs...)

    @info """
    Rotation algorithm converged after $(length(rotation.iterations)) iterations.
          algorithm: $(typeof(method))
          criterion: $(last(rotation.iterations).Q)
    """

    return rotation.L
end

"""
    rotate(Λ, method::RotationMethod; kwargs...)

Perform a rotation of the factor loading matrix Λ and overwrite Λ with the rotated loading
matrix.

See also [`rotate`](@ref).
"""
function rotate!(Λ, method; kwargs...)
    Λ .= rotate(Λ, method; kwargs...)
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
mutable struct RotationState{T<:AbstractMatrix}
    init::T
    A::T
    T::T
    L::T
    iterations::Vector{IterationState}
end

RotationState(init, A) = RotationState(init, A, init, A * init, IterationState[])

"""
    _rotate(A::AbstractMatrix, method::RotationMethod{Orthogonal}; kwargs...)

Implements the algorithm for factor rotation described in Bernaard & Jennrich (2005) for
orthogonal factor rotation.
"""
function _rotate(
    A::AbstractMatrix{TV},
    method::RotationMethod{Orthogonal};
    atol = 1e-6,
    alpha = 1,
    maxiter1 = 1000,
    maxiter2 = 10,
    init::Union{Nothing,AbstractMatrix} = nothing,
) where {TV<:Real}
    state = initialize(init, A)
    Q, ∇Q = criterion_and_gradient(method, state.L)

    # setup intermediate values to avoid unnecessary allocations
    ft = Q
    G = A' * ∇Q
    Tt = similar(state.T)
    Gp = similar(state.T)

    s = zero(eltype(G))

    for i in 1:maxiter1
        project!(Gp, state.T, G)
        s = norm(Gp)

        isconverged(s, atol) && break

        alpha *= 2

        for _ in 1:maxiter2
            X = state.T - alpha * Gp
            @unpack U, Vt = svd(X)
            Tt = U * Vt
            mul!(state.L, A, Tt)
            Q, ∇Q = criterion_and_gradient(method, state.L)
            if (Q < ft - 0.5 * s^2 * alpha)
                state.T = Tt
                ft = Q
                mul!(G, A', ∇Q)
                break
            else
                alpha /= 2
            end
        end

        iteration_state = IterationState(alpha, maxiter2, Q)
        push!(state.iterations, iteration_state)
    end

    isconverged(s, atol) || error("not converged")

    return state
end

"""
    initialize(init, A)

Initialize a [`RotationState`](@ref) with initial values `init` and original loading matrix
`A`. If `init = nothing`, the identity matrix will be used as initial values.
"""
function initialize(init, A::AbstractMatrix{TV}) where {TV}
    p, k = size(A)

    if isnothing(init)
        T = Matrix{TV}(I, k, k)
    else
        T = init
    end

    if size(T) != (k, k)
        throw(ArgumentError("matrix of starting values must be of size ($k, $k)"))
    end

    return RotationState(T, A)
end

"""
    project!(Gp, T, G)

Compute the projection `Gp` of `G` for orthogonal matrices.
"""
function project!(Gp, T, G)
    M = T' * G
    S = (M + M') ./ 2
    Gp .= G - T * S
    return Gp
end

"""
    isconverged(s, atol)

determines the convergence status of the algorithm.
"""
isconverged(s, atol) = s < atol
