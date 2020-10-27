"""
AbstractBlade.jl defines the AbstractBlade type and basic functions

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

export AbstractBlade

"""
    AbstractBlade{<:AbstractFloat}

Supertype for all blade types.

For the AbstractBlade type, the norm and orientation are encoded by the `volume`
of the blade.

Interface
---------

Note: the return value of all methods should preserve the precision of the
AbstractBlade instance (when possible).

### Methods

    grade(B::AbstractBlade)::Int
    basis(B::AbstractBlade; normalized::Bool=true)::Matrix{AbstractFloat}
    volume(B::AbstractBlade{T)::T where {T<:AbstractFloat}
    sign(B::AbstractBlade)::Int8

### Unary Operators

    reciprocal(B::AbstractBlade)::AbstractBlade

### Binary Operators

    proj(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
"""
abstract type AbstractBlade{T<:AbstractFloat} <: AbstractMultivector{T} end

# --- AbstractMultivector interface functions for AbstractBlade type

import Base.sign
import LinearAlgebra.norm

export grades, blades

# Notes: dim() are implemented by subtypes of AbstractBlade.

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
    getindex(B::AbstractBlade, k::Int)

Return the `k`-vector component of blade `B`. When grade(`B`) is equal to `k`,
return a Vector containing `B`. Otherwise, return an empty vector.
"""
Base.getindex(B::AbstractBlade, k::Int) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()

"""
    norm(B::AbstractBlade)

Return the norm of `B`.
"""
norm(B::AbstractBlade) = abs(volume(B))

# --- AbstractBlade interface functions

"""
    sign(B::AbstractBlade)::Int8

Return the sign of `B` relative to its unit basis.
"""
sign(B::AbstractBlade)::Int8 = sign(volume(B))
