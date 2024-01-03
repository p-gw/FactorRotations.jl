"""
    RotationType

An abstract type representing a type of factor rotation.
"""
abstract type RotationType end

"""
    Orthogonal

A type representing an orthogonal rotation type.
"""
struct Orthogonal <: RotationType end

"""
    Oblique

A type representing an oblique rotation type.
"""
struct Oblique <: RotationType end
