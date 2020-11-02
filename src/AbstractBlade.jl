"""
AbstractBlade.jl defines the AbstractBlade type, interface, and core methods

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

export AbstractBlade

# ------ Functions

# Attributes
import Base.sign
export basis, grade, volume

# Unary Operations
export reciprocal

# Binary Operations
export reject

# --- Type definitions

"""
    AbstractBlade{<:AbstractFloat}

Supertype for all blade types.

For the AbstractBlade type, the norm and orientation are encoded by the `volume`
of the blade.

Interface
=========

Note: the return value of all methods should preserve the precision of its
AbstractBlade arguments (when possible).

Functions
---------

### Attributes

    grade(B::AbstractBlade)::Int

    basis(B::AbstractBlade; normalized::Bool=true)::Matrix{AbstractFloat}

    volume(B::AbstractBlade{T)::T where {T<:AbstractFloat}

    sign(B::AbstractBlade)::Int8

### Operations

    dual(B::AbstractBlade, C::AbstractBlade)::AbstractBlade

    reciprocal(B::AbstractBlade)::AbstractBlade

    reject(vectors::Matrix, B::AbstractBlade, normalize::Bool=false)::Matrix
"""
abstract type AbstractBlade{T<:AbstractFloat} <: AbstractMultivector{T} end

# --- Method definitions
#
# Note: the following method definitions are no-op place holders to provide
#       a central location for docstrings.
#

"""
    grade(B)

TODO
"""
grade(B::AbstractBlade) = nothing

"""
    basis(B)

TODO
"""
basis(B::AbstractBlade) = nothing

"""
    volume(B)

TODO
"""
volume(B::AbstractBlade) = nothing

"""
    sign(B::AbstractBlade)::Int8

Return the sign of `B` relative to its unit basis.
"""
sign(B::AbstractBlade)::Int8 = sign(volume(B))

"""
    dual(B)

TODO
"""
dual(B::AbstractBlade) = nothing

"""
    reciprocal(B)

Compute the multiplicative inverse of blade `B`.
"""
reciprocal(B::AbstractBlade) = nothing

# --- Method definitions for AbstractMultivector interface functions

"""
    grades(B::AbstractBlade)

Return a Vector containing the grade of the blade `B`.
"""
grades(B::AbstractBlade) = [grade(B)]

"""
    blades(B::AbstractBlade)

Return a Vector containing the blade `B`.
"""
blades(B::AbstractBlade) = Vector{AbstractBlade}([B])

"""
    norm(B::AbstractBlade)

Return the norm of `B`.
"""
norm(B::AbstractBlade) = abs(volume(B))

"""
    getindex(B::AbstractBlade, k::Int)

Return the `k`-vector component of blade `B`. When grade(`B`) is equal to `k`,
return a Vector containing `B`. Otherwise, return an empty vector.
"""
Base.getindex(B::AbstractBlade, k::Int) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()
