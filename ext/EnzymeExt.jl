module EnzymeExt

using FactorRotations
using Enzyme: gradient!, Reverse

# wrapper for Enzyme.gradient!
FactorRotations.autodiff_gradient!(_::FactorRotations.AutodiffBackend{:Enzyme}, ∇Q::AbstractMatrix, method::RotationMethod, Λ::AbstractMatrix) =
    gradient!(Reverse, ∇Q, Base.Fix1(criterion, method), Λ)

end
