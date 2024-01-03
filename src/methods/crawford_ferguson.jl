"""
    CrawfordFerguson(; kappa, orthogonal = false)

The family of Crawford-Ferguson rotation methods.

## Keyword arguments
- `kappa`: T

## Details
The Crawford-Ferguson family allows both orthogonal and oblique rotation of the
`p`-by`k`-factor loading matrix. If orthogonal rotation is performed, Crawford-Ferguson is
equivalent to Oblimin rotation given the following values for `kappa`:

- `kappa = 0`: Quartimax
- `kappa = 1/p`: Varimax
- `kappa = k/2p`: Equamax
- `kappa = (k - 1)/(p + k - 2)`: Parsimax
- `kappa = 1`: Factor parsimony

"""
struct CrawfordFerguson{T,V} <: RotationMethod{T}
    κ::V
    function CrawfordFerguson(; kappa, orthogonal = false)
        0 <= kappa <= 1 || throw(ArgumentError("kappa must be between 0 and 1"))
        T = orthogonal ? Orthogonal : Oblique
        return new{T,typeof(kappa)}(kappa)
    end
end

function criterion_and_gradient(method::CrawfordFerguson, Λ::AbstractMatrix{T}) where {T}
    @unpack κ = method
    p, k = size(Λ)

    N = ones(T, k, k) |> zerodiag!
    M = ones(T, p, p) |> zerodiag!

    Λsq = Λ .^ 2

    ΛsqN = Λsq * N
    MΛsq = M * Λsq

    Q = (1 - κ) * tr(Λsq' * ΛsqN) / 4 + κ * tr(Λsq' * MΛsq) / 4
    ∇Q = (1 - κ) * Λ .* ΛsqN + κ * Λ .* MΛsq
    return Q, ∇Q
end
