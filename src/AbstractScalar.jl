"""
AbstractScalar.jl defines the AbstractScalar type, interface, and core methods

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports types and interface methods

# ------ Types

export AbstractScalar

# ------ Functions

# Attributes
export value

# --- Type definitions

"""
    AbstractScalar{<:AbstractFloat}

Supertype for all scalar types.

Interface
=========

Note: the return value of all methods should preserve the precision of its
AbstractScalar arguments (when possible).

Functions
---------

### Attributes

    value(B::AbstractScalar{T})::T where {T<:AbstractFloat}
"""
abstract type AbstractScalar{T<:AbstractFloat} <: AbstractBlade{T} end

# --- Method definitions for AbstractBlade interface functions

"""
    grade(B::AbstractScalar)

Return 0.
"""
grade(B::AbstractScalar) = 0

"""
    basis(B::AbstractScalar; normalized::Bool=true)

Return 1.
"""
basis(B::AbstractScalar; normalized::Bool=true) =
    normalized ? 1 : value(B)

"""
    volume(B::AbstractScalar)

Return the value of `B`.
"""
volume(B::AbstractScalar) = value(B)

-(B::AbstractScalar) = Scalar(-value(B))

reciprocal(B::AbstractScalar) = 1 / B

# --- Method definitions for AbstractMultivector interface functions

"""
    dim(B::AbstractScalar)

Return 0.
"""
dim(B::AbstractScalar) = 0

"""
    reverse(B::AbstractScalar)

Return `B` (the reverse of a scalar is itself).
"""
Base.reverse(B::AbstractScalar) = B

"""
    dual(B::AbstractScalar; dim::Integer)

Compute the dual of `B`. Note that the dimension of the embedding space must
be explicitly specified.
"""
function dual(B::AbstractScalar; dim::Union{Integer, Nothing}=nothing)
    if isnothing(dim)
        error("The dual of a scalar is not well-defined if `dim` is not " *
              "specified")
    end

    mod(dim, 4) < 2 ?
        Pseudoscalar(dim, value(B)) :
        Pseudoscalar(dim, -value(B))
end

# --- Comparison methods

# ------ ==(B, C)

# B::AbstractScalar, C::AbstractScalar
==(B::AbstractScalar, C::AbstractScalar) = (value(B) == value(C))

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (value(B) == x)

# ------ isapprox(B, C)

# B::AbstractScalar, C::AbstractScalar
isapprox(B::AbstractScalar{T1}, C::AbstractScalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    isapprox(value(B), value(C), atol=atol, rtol=rtol)

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
isapprox(B::AbstractScalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    isapprox(x, value(B), atol=atol, rtol=rtol)

isapprox(x::Real, B::AbstractScalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    isapprox(B, x, atol=atol, rtol=rtol)

# --- Utility methods

"""
    convert(::Type{S}, B::AbstractScalar) where {T<:AbstractFloat,
                                                 S<:AbstractMultivector{T}}

Convert AbstractScalar to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::AbstractScalar) where {T<:AbstractFloat,
                                             S<:AbstractMultivector{T}} =
    T == typeof(value(B)) ? B : Scalar{T}(value(B))
