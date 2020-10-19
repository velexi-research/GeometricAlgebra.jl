"""
The multivector.jl submodule defines Multivector types.

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
export AbstractMultivector, Multivector

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
abstract type AbstractMultivector end


# Multivector
"""
    struct Multivector{T<:AbstractFloat} <: AbstractMultivector

Multivector represented with the floating-point precision of type `T`.
"""
struct Multivector{T<:AbstractFloat} <: AbstractMultivector
    #=
      Fields
      ------
      * `summands`: collection of blades that sum to the value of the
        multivector

      * `norm`: norm of multivector

      * `reduced`: True if multivector is guaranteed to be reduced to a sum
        of orthogonal blades; False otherwise.
    =#
    summands::SortedDict{Int, Vector{AbstractBlade}}
    norm::T
    reduced::Bool

    """
        Multivector{T}(blades::Vector{AbstractBlade};
                       reduced::Bool=false) where {T<:AbstractFloat}

    Construct a Multivector from a of vector of blades. When `reduced` is true,
    the summands are combined so that the multivectors of grade ``k`` form
    an orthogonal basis for the subspace of ``k``-vectors.
    """
    function Multivector{T}(blades::Vector{<:AbstractBlade};
                            reduced::Bool=false) where {T<:AbstractFloat}

        # --- Construct Multivector

        # Construct `summands`. Sort blades by grade.
        summands = SortedDict{Int, Vector{AbstractBlade}}()
        for B in blades
            k = grade(B)
            k_vectors = get!(summands, k, [])
            push!(k_vectors, convert(AbstractBlade{T}, B))
        end

        # Reduce multivector
        if reduced
            # TODO: reduce k-vectors
        end

        # Compute norm
        norm = 0  # TODO: implement

        # Return new Multivector
        new(summands, norm, false)
    end
end

"""
    Multivector(blades::Vector{<:AbstractBlade}; reduced::Bool=false)

Construct a Multivector from a of vector of blades. When `reduced` is true, the
summands are combined so that the multivectors of grade ``k`` form an
orthogonal basis for the subspace of ``k``-vectors.

The precision of the Multivector is inferred from the precision of the first
element of `blades`.
"""
Multivector(blades::Vector{<:AbstractBlade}; reduced::Bool=false) =
    Multivector{typeof(volume(blades[1]))}(blades, reduced=reduced)


# --- Basic Multivector functions

# Exports
export grades, summands

"""
    grades(M::Multivector; collect=true)

Return the grades of the blades that the multivector is composed of. When
`collect` is `false`, an iterator over the keys is returned.
"""
grades(M::Multivector; collect=true) =
    collect ? Base.collect(keys(M.summands)) : keys(M.summands)

"""
    getindex(M::Multivector, grade::Integer)

Return multivectors of `M` with the specified `grade`.
"""
getindex(M::Multivector, grade::Integer) =
    grade in grades(M) ?
        M.summands[grade] :
        error("`M` does not contain any multivectors of grade $grade")

"""
    summands(M::Multivector)

Return the blades that the multivector is composed of.
"""
summands(M::Multivector) = M.summands

"""
    norm(M::Multivector)

Return the norm of the multivector.
"""
norm(M::Multivector) = M.norm


# --- Utility functions

# Imports
import Base.convert

"""
    convert(::Type{S}, M::Multivector) where {T<:AbstractFloat,
                                              S<:Multivector{T}}

Convert Blade to have the floating-point precision of type `T`.
"""
convert(::Type{S}, M::Multivector) where {T<:AbstractFloat, S<:Multivector{T}} =
    T == typeof(norm(M)) ?
        M : Multivector{T}(reduce(vcat, values(summands(M))))
