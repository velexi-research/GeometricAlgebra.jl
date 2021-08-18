"""
AbstractMultivector.jl defines the AbstractMultivector type and interface

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

# ------ Types

export AbstractMultivector

# ------ Functions

# Attributes
import Base.getindex
import LinearAlgebra.norm
export dim, grades, blades, norm

# Unary operations
import Base.:(-), Base.reverse
export inverse, dual

# Comparison operations
import Base.:(==), Base.isapprox
import Base.iszero, Base.isone

# Utility functions
import Base.zero, Base.one
import Base.convert

# --- Type definitions

"""
    AbstractMultivector{<:AbstractFloat}

Supertype for all multivector types.

Interface
=========

Attributes
----------

    dim(M::AbstractMultivector)::Int

    grades(M::AbstractMultivector)::Vector{Int}

    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}

    norm(M::AbstractMultivector)::AbstractFloat

    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

Unary Operations
----------------

    inverse(M::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector)::AbstractMultivector

    reverse(M::AbstractMultivector)::AbstractMultivector

    dual(M::AbstractMultivector)::AbstractMultivector

Binary Operations
-----------------

    +(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    -(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    *(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    /(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    wedge(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    ∧(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    contract_left(M::AbstractMultivector,
                  N::AbstractMultivector)::AbstractMultivector
    <(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    contract_right(M::AbstractMultivector,
                   N::AbstractMultivector)::AbstractMultivector
    >(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    dot(M::AbstractMultivector, N::AbstractMultivector;
        left=true)::AbstractMultivector
    ⋅(M::AbstractMultivector, N::AbstractMultivector;
      left=true)::AbstractMultivector

    project(M::AbstractMultivector, B::AbstractBlade)::AbstractMultivector

Comparison Functions
--------------------

    ==(M::AbstractMultivector, N::AbstractMultivector)::Bool
    isequal(M::AbstractMultivector, N::AbstractMultivector)::Bool

    ≈(M::AbstractMultivector, N::AbstractMultivector)::Bool
    isapprox(M::AbstractMultivector, N::AbstractMultivector)::Bool

Utility Functions
-----------------

    convert(::Type{T}, M::AbstractMultivector) where {T<:AbstractFloat}

Implementation
==============

* The return value of all methods should preserve the precision of its
  AbstractMultivector arguments (when possible).

"""
abstract type AbstractMultivector{T<:AbstractFloat} end

# --- Method definitions
#
# Note: the following method definitions are no-op place holders and intended
#       to be extended.

"""
    dim(M::AbstractMultivector)::Int

Return dimension of the real space that `M` is embedded within.
"""
function dim end

"""
    grades(M::Multivector; collect=true)::Vector{Int}

Return the grades of the nonzero `k`-vector components of `M`. When `collect`
is `false`, an iterator over the grades is returned.
"""
function grades end

"""
    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}

Return the blades that `M` is composed of.
"""
function blades end

"""
    norm(M::AbstractMultivector)::AbstractFloat

Compute the norm of `M`.
"""
function norm end

"""
    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

Return the `k`-vector component of `M`.

Notes
=====

* When `M` is a blade, `M[k]` is a Vector containing `M` if the grade of `M`
  is equal to `k`; otherwise, `M[k]` is an empty vector.
"""
function getindex end

"""
    inverse(M::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector)::AbstractMultivector

Compute the additive inverse of `M`.
"""
function inverse end

"""
    inverse(M::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector)::AbstractMultivector

Compute the additive inverse of `M`.
"""
-(M::AbstractMultivector) = inverse(M)

"""
    reverse(M::AbstractMultivector)::AbstractMultivector

Compute the reverse of `M`.
"""
function reverse end

"""
    dual(M::AbstractMultivector)::AbstractMultivector

Compute the dual of `M` (relative to the pseudoscalar of the geometric algebra
that `M` is an element of).
"""
function dual end

# --- Comparison methods

# B::AbstractMultivector, x::Real
# x::Real, B::AbstractMultivector
isapprox(M::AbstractMultivector, x::Real) = false
isapprox(x::Real, M::AbstractMultivector) = false

# B::AbstractMultivector, v::Vector
# v::Vector, B::AbstractMultivector
isapprox(M::AbstractMultivector, v::Vector) = false
isapprox(v::Vector, M::AbstractMultivector) = false

iszero(M::AbstractMultivector) = (M === zero(M))

isone(M::AbstractMultivector) = (M === one(M))

# --- Utility methods

zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector}) = Zero{Float64}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()

one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector}) = One{Float64}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()

"""
    convert(::Type{T}, M::AbstractMultivector) where {T<:AbstractFloat}

Convert AbstractMultivector to have the floating-point precision of type `T`.
"""
function convert end

# --- Non-exported utility functions

"""
    assert_dim_equal(M::AbstractMultivector, N::AbstractMultivector)
    assert_dim_equal(M::AbstractMultivector, v::Vector)
    assert_dim_equal(v::Vector, M::AbstractMultivector)

Assert that the dimensions of `M` and `N` are the same.
"""
function assert_dim_equal(M::AbstractMultivector, N::AbstractMultivector)
    if dim(M) != dim(N)
        throw(DimensionMismatch("`dim(M)` not equal to `dim(N)`"))
    end
end

function assert_dim_equal(M::AbstractMultivector, v::Vector)
    if dim(M) != length(v)
        throw(DimensionMismatch("`dim(M)` not equal to `length(v)`"))
    end
end
assert_dim_equal(v::Vector, M::AbstractMultivector) = assert_dim_equal(M, v)
