struct Cubimax <: RotationMethod{Orthogonal} end

criterion(method::Cubimax, Λ::AbstractMatrix) = -sum(x -> abs(x^3), Λ)
