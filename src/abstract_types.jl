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
    abstract type AbstractMultivector

Supertype for all multivector types.

Methods
-------
    grades(M::AbstractMultivector)::Vector
    summands(M::AbstractMultivector)::SortedDict
    norm(M::AbstractMultivector)::AbstractFloat
    reduce(M::AbstractMultivector)::AbstractMultivector
    getindex(M::AbstractMultivector, grade::Integer)::Vector{<:AbstractBlade}

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
    abstract type AbstractBlade{T<:AbstractFloat}

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

"""
    grades(B::AbstractBlade)

Return a Vector containing the grade of `B`.
"""
grades(B::AbstractBlade) = Vector([grade(B)])

"""
    summands(B::AbstractBlade)

Return the a SortedDict containing the `B`.
"""
summands(B::AbstractBlade) =
    SortedDict{Int, Vector{AbstractBlade}}([(grade(B), [B])])

"""
    getindex(B::AbstractBlade, k::Integer)

When grade(`B`) is equal to `k`, return a Vector containing `B`. Otherwise,
return an empty vector.
"""
getindex(B::AbstractBlade, k::Integer) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()
