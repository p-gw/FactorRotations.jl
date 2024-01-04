"""
    rotate(Λ, method::RotationMethod, kwargs...)

Perform a rotation of the factor loading matrix `Λ`.

## Keyword arguments
- `atol`:
- `alpha`:
- `maxiter1`:
- `maxiter2`:
- `init`: A k-by-k matrix of starting values for the algorithm.
          If `init = nothing` (the default), random starting values will be generated.
"""
function rotate(Λ, alg, args...; kwargs...)
    rotation = _rotate(Λ, alg, args...; kwargs...)

    @info """
    Rotation algorithm converged after $(length(rotation.iterations)) iterations.
          algorithm: $(typeof(alg))
          criterion: $(last(rotation.iterations).Q)
    """

    return rotation.L
end

function rotate!(Λ, args...; kwargs...)
    Λ .= rotate(Λ, args...; kwargs...)
    return Λ
end

struct IterationState{T<:Real,V<:Real}
    alpha::T
    maxiter::Int
    Q::V
end

mutable struct RotationState{T<:AbstractMatrix}
    init::T
    A::T
    T::T
    L::T
    iterations::Vector{IterationState}
end

RotationState(init, A) = RotationState(init, A, init, A * init, IterationState[])

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

"compute the projection `Gp` of `G` for orthogonal matrices"
function project!(Gp, T, G)
    M = T' * G
    S = (M + M') ./ 2
    Gp .= G - T * S
    return Gp
end

isconverged(s, atol) = s < atol
