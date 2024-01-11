function kaiser_normalize!(Λ::AbstractMatrix)
    weights = communalities(Λ)
    Λ ./= weights
    return Λ, weights
end

kaiser_normalize(Λ) = kaiser_normalize!(copy(Λ))

function kaiser_denormalize!(Λ::AbstractMatrix, weights::AbstractVector)
    Λ .*= weights
    return Λ
end

kaiser_denormalize(Λ, weights) = kaiser_denormalize!(copy(Λ), weights)

communalities(Λ) = norm.(eachrow(Λ))
