"""
    FactorRotation{T <: Real}

A type holding results of a factor rotation.

## Fields
- `L`: The rotated factor loading matrix
- `T`: The factor rotation matrix
- `phi`: The factor correlation matrix
- `weights`: Normalization weights
"""
struct FactorRotation{T}
    L::Matrix{T}
    T::Matrix{T}
    phi::Matrix{T}
    weights::Vector{T}
end

function FactorRotation(L, T, weights)
    return FactorRotation(L, T, T' * T, weights)
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
- `alpha`: Sets the inital value for alpha (default: 1).
- `f_atol`: Sets the absolute tolerance for the comparison of minimum criterion values when
            with random starts (default: 1e-6).
- `g_atol`: Sets the absolute tolerance for convergence of the algorithm (default: 1e-6).
- `init`: A k-by-k matrix of starting values for the algorithm.
          If `init = nothing` (the default), the identity matrix will be used as starting
          values.
- `maxiter1`: Controls the number of maximum iterations in the outer loop of the algorithm
              (default: 1000).
- `maxiter2`: Controls the number of maximum iterations in the inner loop of the algorithm
              (default: 10).
- `normalize`: Perform Kaiser normalization before rotation of the loading matrix
               (default: false).
- `randomstarts`: Determines if the algorithm should be started from random starting values.
                  If `randomstarts = false` (the default), the algorithm is calculated once
                  for the initial values provided by `init`.
                  If `randomstarts = true`, the algorithm is started 100 times from random
                  starting matrices.
                  If `randomstarts = x::Int`, the algorithm is started `x` times from random
                  starting matrices.
- `reflect`: Switch signs of the columns of the rotated loading matrix such that the sum of
             loadings is non-negative for all columns (default: true)
- `use_threads`: Parallelize random starts using threads (default: false)
- `verbose`: Print logging statements (default: true)
- `logperiod`: How frequently to report the optimization state (default: 100).

## Return type
The `rotate` function returns a [`FactorRotation`](@ref) object.
If `randomstarts` were requested, then `rotate` returns the [`FactorRotation`](@ref) object
with minimum criterion value.

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
function rotate(
    Λ,
    method;
    verbose = VERBOSITY[],
    randomstarts = false,
    normalize = false,
    reflect = true,
    f_atol = 1e-6,
    g_atol = 1e-6,
    use_threads = false,
    kwargs...,
)
    loglevel = verbose ? Logging.Info : Logging.Debug
    starts = parse_randomstarts(randomstarts)

    # pre-processing
    if normalize
        @logmsg loglevel "Performing Kaiser normalization of loading matrix."
        L, weights = kaiser_normalize(Λ)
    else
        L, weights = Λ, ones(eltype(Λ), size(Λ, 1))
    end

    # rotation
    if starts == 0
        rotation = _rotate(L, method; g_atol, loglevel, kwargs...)
    else
        if :init in keys(kwargs)
            @warn "Requested random starts but keyword argument `init` was provided. Ignoring initial starting values in `init`."
        end

        f = use_threads ? Folds.map : map
        states = f(1:starts) do _
            init = random_orthogonal_matrix(size(L, 2))
            return _rotate(L, method; g_atol, loglevel, kwargs..., init)
        end

        rotation = argmin(minimumQ, states)
        Q_mins = minimumQ.(states)
        Q_min = minimumQ(rotation)
        n_at_Q_min = sum(isapprox(Q, Q_min) for Q in Q_mins)
        n_diverged = sum(is_diverged(s) for s in states)

        @logmsg loglevel "Finished $(starts) rotations with random starts."

        if n_diverged == starts
            @warn "All $(starts) rotations did not converge. Please check the provided rotation method and/or loading matrix."
        elseif n_diverged > 0
            @warn "There were $(n_diverged) rotations that did not converge. Please check the provided rotation method and/or loading matrix."
        else
            @logmsg loglevel "There were 0 rotations that did not converge."
        end
        @logmsg loglevel "$(n_at_Q_min) rotations converged to the same minimum value, Q = $(Q_min)"
    end

    # post-processing
    if normalize
        @logmsg loglevel "Denormalizing rotated loading matrix."
        kaiser_denormalize!(rotation.L, weights)
    end

    rot = FactorRotation(rotation.L, rotation.T, weights)

    reflect && reflect!(rot)

    return rot
end

function parse_randomstarts(x::Bool; default = 100)
    starts = x ? default : 0
    return starts
end

function parse_randomstarts(x::Int)
    x > 0 && return x
    msg = "Invalid value argument $(x) for `randomstarts`. Please provide an integer > 0 or set `randomstarts = true` to use the default."
    throw(ArgumentError(msg))
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
mutable struct RotationState{RT<:RotationType,T<:AbstractMatrix,S<:IterationState}
    init::T
    A::T
    T::T
    Ti::Union{Nothing,T}
    L::T
    iterations::Vector{S}
    is_converged::Bool
end

function RotationState(::Type{T}, init, A) where {T<:RotationType}
    if T <: Orthogonal
        Ti = nothing
        L = A * init
    elseif T <: Oblique
        Ti = inv(init)
        L = A * Ti'
    else
        throw(ArgumentError("Unsupported rotation type $(T)"))
    end
    S = IterationState{eltype(A),eltype(A)}
    return RotationState{T,typeof(A),S}(init, A, init, Ti, L, S[], false)
end

is_converged(state::RotationState) = state.is_converged
is_diverged(state::RotationState) = !is_converged(state)

function minimumQ(state::RotationState)
    return length(state.iterations) > 0 ? last(state.iterations).Q : NaN
end

"""
    _rotate(A::AbstractMatrix, method::RotationMethod{Orthogonal}; kwargs...)

Implements the algorithm for factor rotation described in Bernaard & Jennrich (2005) for
orthogonal factor rotation.
"""
function _rotate(
    A::AbstractMatrix{TV},
    method::RotationMethod{RT};
    g_atol = 1e-6,
    alpha = 1,
    maxiter1 = 1000,
    maxiter2 = 10,
    init::Union{Nothing,AbstractMatrix} = nothing,
    logperiod::Integer = 100,
    loglevel,
) where {RT,TV<:Real}
    @logmsg loglevel "Initializing rotation using algorithm $(typeof(method))."
    state = initialize(RT, init, A; loglevel)
    ∇Q = similar(state.L)
    Q = criterion_and_gradient!(∇Q, method, state.L)

    @logmsg loglevel "Initial criterion value = $(Q)"

    # preallocate variables to avoid unnecessary allocations
    ft = Q
    G = gradient_f!(similar(state.T), state, ∇Q)
    X = similar(state.T)
    Tt = similar(state.T)
    Gp = similar(state.T)
    s = zero(eltype(G))

    α = Float64(alpha) # for type stability

    @logmsg loglevel "Starting optimization..."
    for i in 1:maxiter1
        project_G!(Gp, state, G)
        s = norm(Gp)

        is_converged(s, g_atol) && break

        α *= 2

        for _ in 1:maxiter2
            copy!(X, state.T)
            axpy!(-α, Gp, X)
            project_X!(Tt, state, X)
            update_state!(state, Tt)

            Q = criterion_and_gradient!(∇Q, method, state.L)

            if (Q < ft - 0.5 * s^2 * α)
                # update state.T (and reuse the old one for the next iteration)
                Tt, state.T = state.T, Tt
                ft = Q
                gradient_f!(G, state, ∇Q)
                break
            else
                α /= 2
            end
        end

        if (i == 1 || i == maxiter1 || mod(i, logperiod) == 0)
            @logmsg loglevel "Current optimization state:" iteration = i criterion = Q alpha =
                α
        end

        iteration_state = IterationState(α, maxiter2, Q)
        push!(state.iterations, iteration_state)
    end

    if !is_converged(s, g_atol)
        @warn "Algorithm did not converge after $(maxiter1) iterations (|∇G|=$(s) > $(g_atol))"
    else
        state.is_converged = true
        @logmsg loglevel "Rotation algorithm converged after $(length(state.iterations)) iterations."
        @logmsg loglevel "Final criterion value = $(ft)"
    end

    return state
end

"""
    initialize(init, A)

Initialize a [`RotationState`](@ref) with initial values `init` and original loading matrix
`A`. If `init = nothing`, the identity matrix will be used as initial values.
"""
function initialize(
    ::Type{RT},
    init,
    A::AbstractMatrix{TV};
    loglevel,
) where {RT<:RotationType,TV}
    _, k = size(A)

    if isnothing(init)
        @logmsg loglevel "No initial values provided. Using identity matrix as starting value."
        T = Matrix{TV}(I, k, k)
    else
        T = init
    end

    if size(T) != (k, k)
        throw(ArgumentError("matrix of starting values must be of size ($k, $k)"))
    end

    return RotationState(RT, T, A)
end

function gradient_f!(G::AbstractMatrix, state::RotationState{Orthogonal}, ∇Q)
    @unpack A = state
    mul!(G, A', ∇Q)
    return G
end

function gradient_f!(G::AbstractMatrix, state::RotationState{Oblique}, ∇Q)
    @unpack L, Ti = state
    mul!(G, Ti', ∇Q' * L, -1, 0)
    return G
end

"""
    project!(Gp, T, G)

Compute the projection `Gp` of `G` and store the results in `Gp`.
"""
function project_G!(Gp, state::RotationState{Orthogonal}, G)
    @unpack T = state
    M = T' * G
    S = M + M'
    mul!(Gp, T, S, -0.5, 0)
    axpy!(1, G, Gp)
    return Gp
end

function project_G!(Gp, state::RotationState{Oblique}, G)
    @unpack T = state
    TG = T .* G
    mul!(Gp, T, diagm(vec(sum(TG, dims = 1))), -1, 0)
    axpy!(1, G, Gp)
    return Gp
end

"""
compute the projection `Tt` of `X`.
"""
function project_X!(Tt, state::RotationState{Orthogonal}, X)
    @unpack U, Vt = svd(X)
    mul!(Tt, U, Vt)
    return Tt
end

function project_X!(Tt, state::RotationState{Oblique}, X)
    v = inv.(sqrt.(sum(abs2, X, dims = 1)))
    mul!(Tt, X, diagm(vec(v)))
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
    is_converged(s, g_atol)

determines the convergence status of the algorithm.
"""
is_converged(s, g_atol) = s < g_atol
