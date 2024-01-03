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
function rotate(A, args...; kwargs...)
    L = copy(A)
    return rotate!(L, args...; kwargs...)
end

function rotate!(
    A::AbstractMatrix{TV},
    method::RotationMethod{Orthogonal};
    atol = 1e-6,
    alpha = 1,
    maxiter1 = 1000,
    maxiter2 = 10,
    init::Union{Nothing,AbstractMatrix} = nothing,
) where {TV<:Real}
    p, k = size(A)

    if isnothing(init)
        T = Matrix(qr(rand(TV, k, k)).Q)
    else
        size(init) == (k, k) ||
            throw(ArgumentError("matrix of starting values must be of size ($k, $k)"))
        T = init
    end

    L = A * T
    Q, ∇Q = criterion_and_gradient(method, L)
    ft = Q
    G = A' * ∇Q

    s = zero(eltype(G))
    maxiter = 0

    for i in 1:maxiter1
        M = T' * G
        S = (M + M') / 2
        Gp = G - T * S

        s = norm(Gp)

        alpha *= 2

        for _ in 1:maxiter2
            X = T - alpha * Gp
            @unpack U, Vt = svd(X)
            Tt = U * Vt
            L = A * Tt
            Q, ∇Q = criterion_and_gradient(method, L)
            if (Q < ft - 0.5 * s^2 * alpha)
                T = Tt
                ft = Q
                G = A' * ∇Q
                break
            else
                alpha /= 2
            end
        end

        maxiter = i
        s < atol && break
    end

    s < atol || error("not converged")

    @info """
    Rotation algorithm converged after $maxiter iterations.
          Algorithm: $(method)
          criterion: $(ft)
    """

    A .= L
    return A
end
