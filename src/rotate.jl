"""
    rotate(Λ, method::RotationMethod, kwargs...)

Perform a rotation of the factor loading matrix `Λ`.

## Keyword arguments
- `atol`:
- `alpha`:
- `maxiter1`:
- `maxiter2`:
- `init`: A p-by-p matrix of starting values for the algorithm.
          If `init = nothing` (the default), random starting values will be generated.
"""
rotate(args...; kwargs...) = _rotate(args...; kwargs...)

# https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=f448431f1d1e2179cf1db64b86b03cc9bbe2ab2f
function _rotate(
    A::AbstractMatrix{TV},
    method::RotationMethod{Orthogonal};
    atol = 1e-6,
    alpha = 1,
    maxiter1 = 1000,
    maxiter2 = 10,
    init = nothing,
) where {TV<:Real}
    if isnothing(init)
        p = size(A, 2)
        T = qr(rand(TV, p, p)).Q
    else
        T = init
    end
    T = Matrix{Float64}(I, size(A, 2), size(A, 2))

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

    return L
end
