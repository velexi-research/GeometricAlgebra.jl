"""
Pseudoscalar.jl defines the Pseudoscalar type and basic functions

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

export Pseudoscalar

"""
    struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade

Pseudoscalar ((n-1)-blade) represented with the floating-point precision of
type `T`. The `basis` for a Pseudoscalar is the standard basis for ``R^n``.
The norm and orientation of a Pseudoscalar are encoded in its `value`. The
norm of a Pseudoscalar is equal to `abs(value)` and the orientation of a
Pseudoscalar relative to the standard basis for ``R^n`` is equal to
`sign(value)`.
"""
struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade{T}
    #=
      Fields
      ------
      * `dim`: the dimension of the space that the blade is embedded in

      * `value`: the value of the pseudoscalar
    =#
    dim::Int
    value::T

    """
    Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having
    the specified `value`.
    """
    Pseudoscalar{T}(dim::Integer, value::Real) where {T<:AbstractFloat} =
        dim <= 0 ? error("`dim` must be positive") :
            value == 0 ? zero(Pseudoscalar{T}) : new(dim, T(value))
end

"""
    Pseudoscalar(dim::Integer, value::Real)

Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having the
specified `value`.

When the precision is not specified, the following rules are applied to set
the precision of the Pseudoscalar.

* If `value` is a floating-point value, the precision of the Pseudoscalar is
  inferred from the precision of `value`.

* If `value` is an integer, the precision of the Pseudoscalar defaults to
  `Float64`.
"""
Pseudoscalar(dim::Integer, value::AbstractFloat) =
    Pseudoscalar{typeof(value)}(dim, value)

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))


"""
    Pseudoscalar(B::Pseudoscalar{T};
                 value::Real=value(B)) where {T<:AbstractFloat}

Copy constructor. Construct a Pseudoscalar representing the same space as
`B` having the specified `value`.
"""
Pseudoscalar(B::Pseudoscalar; value::Real=value(B)) =
    Pseudoscalar{typeof(B.value)}(dim(B), value)

# --- AbstractMultivector interface functions for Pseudoscalar type

export dim

"""
    dim(B::Pseudoscalar)

Return dimension of space that `B` is embedded in.
"""
dim(B::Pseudoscalar) = B.dim

# --- AbstractBlade interface functions for Pseudoscalar type

import LinearAlgebra.I

export grade, basis, volume

"""
    grade(B::Pseudoscalar)

Return the grade of the dimension of the space spanned by `B`.
"""
grade(B::Pseudoscalar) = B.dim

"""
    basis(B::Pseudoscalar)

Return LinearAlgebra.I.
"""
basis(B::Pseudoscalar) = LinearAlgebra.I

"""
    volume(B::Pseudoscalar)

Return the value of `B`.
"""
volume(B::Pseudoscalar) = value(B)

# --- Pseudoscalar functions

export value

"""
    value(B::Pseudoscalar)::Real

Return the value of `B`.
"""
value(B::Pseudoscalar) = B.value

# --- Utility functions

import Base.convert

"""
    convert(::Type{S}, B::Pseudoscalar)
        where {T<:AbstractFloat, S<:Pseudoscalar{T}}

Convert Pseudoscalar to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Pseudoscalar) where {T<:AbstractFloat,
                                           S<:AbstractMultivector{T}} =
    T == typeof(value(B)) ? B : Pseudoscalar{T}(dim(B), value(B))
