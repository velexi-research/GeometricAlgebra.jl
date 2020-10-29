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
# --- Exports

# Types
export Multivector

# --- Imports

using DataStructures: SortedDict

"""
    struct Multivector{T<:AbstractFloat} <: AbstractMultivector

Multivector represented with the floating-point precision of type `T`.
"""
struct Multivector{T<:AbstractFloat} <: AbstractMultivector{T}
    #=
      Fields
      ------
      * `parts`: collection of k-vectors that sum to the multivector

      * `norm`: norm of the multivector
    =#
    dim::Int
    parts::SortedDict{Int, Vector{AbstractBlade}}
    norm::T

    """
    Construct a Multivector from a of vector of blades. For each grade ``k``,
    the blades used to represent the ``k``-vector part of the multivector form
    an orthogonal basis for the subspace of ``k``-vectors.
    """
    function Multivector{T}(
            blades::Vector{<:AbstractBlade}) where {T<:AbstractFloat}

        # --- Handle edge cases

        # Return Scalar(0) if number of blades is 0.
        if length(blades) == 0
            return zero(Blade{T})
        end

        # --- Construct Multivector

        dim_ = dim(blades[1])

        # Construct `parts`. Sort blades by grade.
        parts = SortedDict{Int, Vector{AbstractBlade}}()
        for B in blades
            # Ignore Scalar(0)
            if volume(B) == 0
                continue
            end

            k_vectors = get!(parts, grade(B), [])
            push!(k_vectors, convert(AbstractBlade{T}, B))
        end

        # Return AbstractBlade if there is only one term
        if length(parts) == 1
            k_vectors = first(values(parts))
            if length(k_vectors) == 1
                return first(k_vectors)
            end
        end

        # --- Reduce multivector

        # Reduce scalar, vector, and pseudoscalar parts
        for k in [0, 1, dim_]
            if k in keys(parts)
                parts[k] = [reduce(+, parts[k])]
            end
        end

        # TODO: reduce remaining k-vectors

        # Compute norm
        norm = 0  # TODO: implement

        # Return new Multivector
        new(dim_, parts, norm)
    end

    """
    Type conversion constructor. Construct a copy of the Multivector converted
    to the specified precision.
    """
    Multivector{T}(M::Multivector) where {T<:AbstractFloat} =
        T == typeof(norm(M)) ?
            M :
            Multivector{T}(reduce(vcat, blades(M)))
end

"""
    Multivector(blades::Vector{<:AbstractBlade})

Construct a Multivector from a of vector of blades. For each grade ``k``, the
blades used to represent the ``k``-vector part of the multivector form an
orthogonal basis for the subspace of ``k``-vectors.

The precision of the Multivector is inferred from the precision of the first
element of `blades`.
"""
Multivector(blades::Vector{<:AbstractBlade}) =
    length(blades) > 0 ?
        Multivector{typeof(volume(blades[1]))}(blades) :
        Multivector{Float64}(blades)

Multivector(multivectors::Vector{<:AbstractMultivector}) =
    length(multivectors) > 0 ?
        Multivector(reduce(vcat, map(blades, multivectors))) :
        Multivector{Float64}(Vector{AbstractBlade}())



# --- AbstractMultivector interface functions for Multivector type

"""
    grades(M::Multivector; collect=true)

Return the grades of the blades that multivector `M` is composed of. When
`collect` is `false`, an iterator over the grades is returned.
"""
grades(M::Multivector; collect=true) =
    collect ?
        Base.collect(keys(M.parts)) :
        keys(M.parts)

"""
    blades(M::Multivector)

Return the blades that the multivector `M` is composed of.
"""
blades(M::Multivector) = reduce(vcat, values(M.parts))

"""
    getindex(M::Multivector, k::Int)

Return the `k`-vector component of multivector `M`.
"""
Base.getindex(M::Multivector, k::Int) =
    k in grades(M) ?  M.parts[k] : Vector{AbstractBlade}()

"""
    norm(M::Multivector)

Return the norm of the multivector `M`.
"""
norm(M::Multivector) = M.norm

"""
    convert(::Type{S}, M::Multivector) where {T<:AbstractFloat,
                                              S<:Multivector{T}}

Convert Blade to have the floating-point precision of type `T`.
"""
convert(::Type{S}, M::Multivector) where {T<:AbstractFloat, S<:Multivector{T}} =
    Multivector{T}(M)
