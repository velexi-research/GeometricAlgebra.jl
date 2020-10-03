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



# --- Types

# Exports
export AbstractMultivector, Multivector

# AbstractMultivector
"""
    abstract type AbstractMultivector

Supertype for all multivector types.

Methods
-------
    grades(M::AbstractMultivector)::KeySet
    summands(M::AbstractMultivector)::Dict
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

TODO
"""
struct Multivector{T<:AbstractFloat} <: AbstractMultivector
    #=
      Fields
      ------
      * `summands`: collection of blades that sum to the value of the
        multivector

      * `reduced`: True if multivector is guaranteed to be reduced to a sum
        of orthogonal blades; False otherwise.
    =#
    summands::Dict{Int, Vector{AbstractBlade}}
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
        summands = Dict{Int, Vector{AbstractBlade}}()
        for B in blades
            k = grade(B)
            k_vectors = get!(summands, k, [])
            push!(k_vectors, B)
        end

        # Reduce multivector
        if reduced
            # TODO: reduce k-vectors
        end

        # Return new Multivector
        new(summands, false)
    end
end

"""
    Multivector(blades::Vector{T};
          value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

Construct a Multivector from a of vector of blades. When `reduced` is true, the
summands are combined so that the multivectors of grade ``k`` form an
orthogonal basis for the subspace of ``k``-vectors.
"""
Multivector(blades::Vector{<:AbstractBlade}; reduced::Bool=false) =
    Multivector{typeof(volume(blades[1]))}(blades, reduced=reduced)


# --- Basic Multivector functions

# Exports
export grades, summands

"""
TODO
"""
grades(M::Multivector) = keys(M.summands)

"""
TODO
"""
summands(M::Multivector) = M.summands
