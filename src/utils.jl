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

"""
    centercols!(m::AbstractMatrix)

Efficiently substract the column mean for each column in `m`.
"""
function centercols!(m::AbstractMatrix)
    for col in eachcol(m)
        col .-= mean(col)
    end
    return m
end

"""
    ConvergenceError

The algorithm did not converge. `msg` provides a descriptive error message.
"""
struct ConvergenceError <: Exception
    msg::String
end

Base.show(io::IO, err::ConvergenceError) = print(io, "ConvergenceError: $(err.msg)")
