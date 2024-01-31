"""
    Simplimax(; m::Int)

The Simplimax rotation method.
"""
struct Simplimax <: RotationMethod{Oblique}
    m::Int
    function Simplimax(; m::Int)
        m < 1 && throw(ArgumentError("m must be greater or equal to 1"))
        return new(m)
    end
end

function criterion(method::Simplimax, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    λm = nthsmallest(Λsq, method.m)
    Q = tr(Λsq' * (Λsq .<= λm)) / 2
    return Q
end

function criterion_and_gradient(method::Simplimax, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    λm = nthsmallest(Λsq, method.m)
    Λind = Λsq .<= λm

    Q = tr(Λsq' * Λind) / 2
    ∇Q = Λ .* Λind
    return Q, ∇Q
end
