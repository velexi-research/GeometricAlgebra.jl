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
export dim, blades, grades, norm

# Unary operations
import Base.:(-), Base.reverse
export dual

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

Note: the return value of all methods should preserve the precision of its
AbstractMultivector arguments (when possible).

Attributes
----------

    dim(M::AbstractMultivector)::Int

    grades(M::AbstractMultivector)::Vector{Int}

    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}

    norm(M::AbstractMultivector{T})::T where {T<:AbstractFloat}

    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

Unary Operations
----------------

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

    contractl(M::AbstractMultivector,
              N::AbstractMultivector)::AbstractMultivector
    <(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    contractr(M::AbstractMultivector,
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

    isapprox(M::AbstractMultivector, N::AbstractMultivector)::Bool
    ≈(M::AbstractMultivector, N::AbstractMultivector)::Bool

Utility Functions
-----------------

    convert(::Type{T}, M::AbstractMultivector)
        where {T<:AbstractMultivector{<:AbstractFloat}}
"""
abstract type AbstractMultivector{T<:AbstractFloat} end

# --- Method definitions
#
# Note: the following method definitions are no-op place holders and intended
#       to be extended.

"""
    dim(M)

Return dimension of the real space that `M` is embedded within.
"""
function dim end

"""
    blades(M)

Return the blades that `M` is composed of.
"""
function blades end

"""
    grades(M)

Return the grades of the nonzero `k`-vector components of `M`.
"""
function grades end

"""
    norm(M)

Compute the norm of `M`.
"""
function norm end

"""
    getindex(M::AbstractMultivector, k::Int)

Return the `k`-vector component of multivector `M`.
"""
function getindex end

"""
    inverse(M), -(M)

Compute the additive inverse of a multivector `M`.
"""
function inverse end
-(M::AbstractMultivector) = inverse(M)

"""
    reverse(M)

Compute the reverse of a multivector `M`.
"""
function reverse end

"""
    dual(M)

Compute the dual of a multivector `M` (relative to the space that the
geometric algebra is extended from).

TODO: find better choice of words than "extended from"
"""
function dual end

# --- Comparison methods

"""
    ==(M::AbstractMultivector, N::AbstractMultivector)

Return true if `M` and `N` are equal; otherwise, return false.
"""
function == end

"""
    isapprox(M::AbstractMultivector, N::AbstractMultivector)
    ≈(M::AbstractMultivector, N::AbstractMultivector)

Return true if `M` and `N` are approximately equal; otherwise, return false.
"""
function isapprox end

# B::AbstractMultivector, x::Real
# x::Real, B::AbstractMultivector
isapprox(M::AbstractMultivector, x::Real) = false
isapprox(x::Real, M::AbstractMultivector) = false

# B::AbstractMultivector, v::Vector
# v::Vector, B::AbstractMultivector
isapprox(M::AbstractMultivector, v::Vector) = false
isapprox(v::Vector, M::AbstractMultivector) = false

# Comparison to zero()
iszero(M::AbstractMultivector) = (M === zero(M))

# Comparison to one()
isone(M::AbstractMultivector) = (M === one(M))

# --- Utility methods

"""
    zero(M::AbstractMultivector)
    zero(::Type{<:AbstractMultivector})
    zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector}) = Zero{Float64}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()

"""
    one(M::AbstractMultivector)
    one(::Type{<:AbstractMultivector})
    one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector}) = One{Float64}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()

"""
    convert(::Type{T}, M::AbstractMultivector)
        where {T<:AbstractMultivector{<:AbstractFloat}}

Convert AbstractScalar to have the floating-point precision of type `T`.
"""
function convert end

# --- Non-exported utility functions

"""
    assert_dim_equal(M, N)

Assert that the dimensions of `M` and `N` are the same.

Valid arguments
---------------
    assert_dim_equal(M::AbstractMultivector, N::AbstractMultivector)
    assert_dim_equal(M::AbstractMultivector, v::Vector)
    assert_dim_equal(v::Vector, M::AbstractMultivector)
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
