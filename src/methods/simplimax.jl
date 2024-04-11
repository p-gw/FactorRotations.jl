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

function criterion_and_gradient!(∇Q::OptionalGradient, method::Simplimax, Λ::AbstractMatrix)
    Λsq = Λ .^ 2
    λm = nthsmallest(Λsq, method.m)
    Λind = Λsq .<= λm

    Q = dot(Λsq, Λind) / 2
    if !isnothing(∇Q)
        ∇Q .= Λ .* Λind
    end
    return Q
end
