function zerodiag!(x)
    for i in diagind(x)
        x[i] = 0.0
    end
    return x
end

mxlogx(x) = -xlogx(x)

"""
    nthsmallest(m::AbstractArry, n)

Get the nth smallest element of `m`.
"""
nthsmallest(m::AbstractArray, n) = sort(vec(m))[n]

"""
    random_orthogonal_matrix(n::Int)

Return a random orthogonal square matrix of size (n, n).
"""
function random_orthogonal_matrix(n::Int)
    Q, R = qr(randn(n, n))
    O = Q * Diagonal(sign.(diag(R)))
    return O
end
