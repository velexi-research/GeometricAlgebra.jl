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

# Utility functions
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

    proj(M::AbstractMultivector, B::AbstractBlade)::AbstractMultivector

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
dim(M::AbstractMultivector) = nothing

"""
    blades(M)

Return the blades that `M` is composed of.
"""
blades(M::AbstractMultivector) = nothing

"""
    grades(M)

Return the grades of the nonzero `k`-vector components of `M`.
"""
grades(M::AbstractMultivector) = nothing

"""
    norm(M)

Compute the norm of `M`.
"""
norm(M::AbstractMultivector) = nothing

"""
    getindex(M::AbstractMultivector, k::Int)

Return the `k`-vector component of multivector `M`.
"""
Base.getindex(M::AbstractMultivector, k::Int) = nothing

"""
    -(M)

Compute the additive inverse of a multivector `M`.
"""
-(M::AbstractMultivector) = nothing

"""
    reverse(M)

Compute the reverse of a multivector `M`.
"""
reverse(M::AbstractMultivector) = nothing

"""
    dual(M)

Compute the dual of a multivector `M` (relative to the space that the
geometric algebra is extended from).

TODO: find better choice of words than "extended from"
"""
dual(M::AbstractMultivector) = nothing

# --- Comparison methods

"""
    ==(M::AbstractMultivector, N::AbstractMultivector)

Return true if `M` and `N` are equal; otherwise, return false.
"""
==(M::AbstractMultivector, N::AbstractMultivector) = false

"""
    isapprox(M::AbstractMultivector, N::AbstractMultivector)
    ≈(M::AbstractMultivector, N::AbstractMultivector)

Return true if `M` and `N` are approximately equal; otherwise, return false.
"""
isapprox(M::AbstractMultivector, N::AbstractMultivector) = false

# B::AbstractMultivector, x::Real
# x::Real, B::AbstractMultivector
isapprox(M::AbstractMultivector, x::Real) = false
isapprox(x::Real, M::AbstractMultivector) = false

# B::AbstractMultivector, v::Vector
# v::Vector, B::AbstractMultivector
isapprox(M::AbstractMultivector, v::Vector) = false
isapprox(v::Vector, M::AbstractMultivector) = false

# --- Utility methods

"""
    convert(::Type{T}, M::AbstractMultivector)
        where {T<:AbstractMultivector{<:AbstractFloat}}

Convert AbstractScalar to have the floating-point precision of type `T`.
"""
convert(::Type{T}, M::AbstractMultivector) where {S<:AbstractFloat,
                                                  T<:AbstractMultivector{S}} =
    nothing

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
