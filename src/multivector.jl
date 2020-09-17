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
    abstract type AbstractMultivector{T<:AbstractFloat}

Supertype for all multivector types. Multivectors are represented with the
floating-point precision of type `T`.

Methods
-------
    summands(M::AbstractMultivector{T})::Dict
    norm(M::AbstractMultivector{T})::T
    reduce(M::AbstractMultivector{T})::AbstractMultivector{T}

Unary Operations
----------------
    -(M::AbstractMultivector{T})::AbstractMultivector{T}

Binary Operations
------------------
    +(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    *(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    /(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
"""
abstract type AbstractMultivector{T<:AbstractFloat} end


# Multivector
"""
    struct Multivector{T<:AbstractFloat} <: AbstractMultivector{T}

TODO
"""
struct Multivector{T<:AbstractFloat} <: AbstractMultivector{T}
    # Fields
    # ------
    # * `summands`: collection of blades that sum to the value of the
    #   multivector
    #
    # * `reduced`: True if multivector is guaranteed to be reduced to a sum
    #   of orthogonal blades; False otherwise.
    summands::Dict{Int, Vector{AbstractBlade{T}}}
    reduced::Bool

    """
        Multivector{T}(blades::Vector{AbstractBlade{T}};
                       reduced::Bool=false) where {T<:AbstractFloat}

    Construct a Multivector from a of vector of blades. When `reduced` is true,
    the summands are combined so that the multivectors of grade ``k`` form
    an orthogonal basis for the subspace of ``k``-vectors.
    """
    function Multivector{T}(blades::Vector{AbstractBlade{T}};
                            reduced::Bool=false) where {T<:AbstractFloat}

        # --- Construct Multivector

        # Construct `summands`. Sort blades by grade.
        summands = Dict{Int, Vector{AbstractBlade{T}}}()
        for B in blades
            if haskey(summands, grade(B))
                append!(summands[grade(B)], B)
            else
                summands[grade(B)] = Vector{AbstractBlade{T}}([B])
            end
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
    Multivector(vectors::Array{T};
          value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

Construct a Multivector from a of vector of blades. When `reduced` is true, the
summands are combined so that the multivectors of grade ``k`` form an
orthogonal basis for the subspace of ``k``-vectors.
"""
Multivector(blades::Array{T}; reduced::Bool=false) where {T<:AbstractFloat} =
    Multivector{T}(blades, reduced=reduced)

Multivector{T}(blades::Array{<:AbstractFloat};
               reduced::Bool=false) where {T<:AbstractFloat} =
    Multivector{T}(convert.(blades, AbstractBlade{T}), reduced=reduced)


# --- Basic Multivector functions

# Exports
export summands

summands(M::Multivector) = M.summands
