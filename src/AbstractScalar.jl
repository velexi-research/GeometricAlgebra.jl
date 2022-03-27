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

# Types
export AbstractScalar

# Properties
export value

# --- Type definitions

"""
    AbstractScalar{<:AbstractFloat}

Supertype for all scalar types.

Interface
=========

Note: the return value of all methods should preserve the precision of its
AbstractScalar arguments (when possible).

Properties
----------

    value(B::AbstractScalar{T})::T where {T<:AbstractFloat}
"""
abstract type AbstractScalar{T<:AbstractFloat} <: AbstractBlade{T} end

# --- Method definitions
#
# Note: the following method definitions are no-op place holders and intended
#       to be extended.

"""
    value(B::AbstractScalar)::AbstractFloat

Return the value of `B` (with the same precision as `B`).
"""
function value end

# --- Method definitions for AbstractBlade interface functions

grade(B::AbstractScalar) = 0

basis(B::AbstractScalar) = one(value(B))

volume(B::AbstractScalar) = value(B)

inverse(B::AbstractScalar) = Scalar(-value(B))

reciprocal(B::AbstractScalar) = 1 / B

# --- Method definitions for AbstractMultivector interface functions

"""
    dim(B::AbstractScalar)::Int

Return 0.

Notes
=====

* The convention that scalars have zero dimension is adopted because
  (1) scalars exist independently of all geometric algebras and (2) scalars
  are 0-dimensional entities.
"""
dim(B::AbstractScalar) = 0

Base.reverse(B::AbstractScalar) = B

"""
    dual(B::AbstractScalar, dim::Integer)::Pseudoscalar

Compute the dual of `B` relative to a real vector space having dimension `dim`.

Notes
=====

* An error is raised if the dimension of the embedding space is not explicitly
  specified.
"""
dual(B::AbstractScalar, dim::Integer) =
    mod(dim, 4) < 2 ?
        Pseudoscalar(dim, value(B)) :
        Pseudoscalar(dim, -value(B))

dual(B::AbstractScalar) =
    error("The dual of a scalar is not well-defined if `dim` is not " *
          "specified")

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

convert(::Type{T}, B::AbstractScalar) where {T<:AbstractFloat} =
    T == typeof(value(B)) ? B : Scalar{T}(value(B))
