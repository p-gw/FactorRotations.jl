struct Varimax <: RotationMethod{Orthogonal} end

function criterion_and_gradient(::Varimax, Λ::AbstractMatrix{T}) where {T<:Real}
    Λsq = Λ .^ 2
    colmeans = mean(Λsq, dims = 1)
    Λsq_demeaned = Λsq .- colmeans
    Q = -norm(Λsq_demeaned)^2 / 4
    ∇Q = @. -Λ * Λsq_demeaned
    return Q, ∇Q
end
