"""
The abstract_types.jl submodule defines the abstract type hierarchy.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

import DataStructures.SortedDict

# --- Types

# Exports
export AbstractMultivector, AbstractBlade

# AbstractMultivector
"""
    AbstractMultivector{<:AbstractFloat}

Supertype for all multivector types.

Methods
-------
    grades(M::AbstractMultivector)::Vector{Int}
    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}

    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

    norm(M::AbstractMultivector)::AbstractFloat

    TODO: review
    reduce(M::AbstractMultivector)::AbstractMultivector

Unary Operations
----------------
    -(M::AbstractMultivector)::AbstractMultivector

Binary Operations
------------------
    +(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    *(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    /(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
"""
abstract type AbstractMultivector{T<:AbstractFloat} end


# AbstractBlade
"""
    AbstractBlade{<:AbstractFloat}

Supertype for all blade types.

For the AbstractBlade type, the norm and orientation are encoded by the `volume`
of the blade.

Methods
-------
    dim(B)::Int
    grade(B)::Int
    basis(B)::Matrix{AbstractFloat}
    volume(B)::AbstractFloat
    norm(B)::AbstractFloat
    sign(B)::Int8

Unary Operations
----------------
    -(B)::AbstractBlade
    dual(B)::AbstractBlade
    reciprocal(B)::AbstractBlade
    reverse(B)::AbstractBlade

Binary Operations
------------------
    ∧(B, C)::AbstractBlade
    outer(B, C)::AbstractBlade

    ⋅(B, C)::AbstractBlade
    inner(B, C)::AbstractBlade

    *(B, C)::Union{AbstactBlade, AbstractMultivector}
    /(B, C)::Union{AbstactBlade, AbstractMultivector}

    +(B, C)::AbstractMultivector
    -(B, C)::AbstractMultivector

    project(A, B)::AbstractBlade
    dual(A, B)::AbstractBlade
"""
abstract type AbstractBlade{T<:AbstractFloat} <: AbstractMultivector{T} end


# --- AbstractMultivector interface functions: AbstractBlade

export grades, blades

"""
    grades(B::AbstractBlade)

Return a Vector containing the grade of the blade `B`.
"""
grades(B::AbstractBlade) = Vector([grade(B)])

"""
    blades(B::AbstractBlade)

Return a Vector containing the blade `B`.
"""
blades(B::AbstractBlade) = Vector([B])

"""
    getindex(B::AbstractBlade, k::Int)

Return the `k`-vector component of blade `B`. When grade(`B`) is equal to `k`,
return a Vector containing `B`. Otherwise, return an empty vector.
"""
Base.getindex(B::AbstractBlade, k::Int) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()
