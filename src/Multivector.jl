#   Copyright (c) 2020-2022 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Multivector.jl defines the Multivector type and interface
"""

# --- Exports

# ------ Types

export Multivector

# --- Type definitions

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

        # --- Check arguments

        # Check that all non-scalar blades have the same dimension
        dim_ = nothing
        for B in blades
            if !(B isa AbstractScalar)
                if isnothing(dim_)
                    dim_ = dim(B)
                else
                    if dim_ != dim(B)
                        message = "Non-scalar blades in `blades` have " *
                                  "differing dimensions."
                        throw(DimensionMismatch(message))
                    end
                end
            end
        end

        # --- Handle edge cases

        # Return Scalar(0) if number of blades is 0.
        if length(blades) == 0
            return zero(Blade{T})
        end

        # Return a Scalar if all of the blades are scalars
        blades_all_scalar = true
        for B in blades
            if !(B isa AbstractScalar)
                blades_all_scalar = false
                break
            end
        end
        if blades_all_scalar
            return Scalar{T}(value(reduce(+, blades)))
        end

        # --- Construct Multivector

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

        # --- Reduce multivector

        # Reduce scalar, vector, and pseudoscalar parts
        for k in [0, 1, dim_]
            if k in keys(parts)
                parts[k] = [reduce(+, parts[k])]
            end
        end

        # TODO: reduce remaining k-vectors

        # --- Compute norm

        norm_ = 0
        for k_vectors in parts
            norm_ += reduce(+, map(B -> norm(B)^2, k_vectors))
        end
        norm_ = sqrt(norm_)

        # --- Construct result

        # Return AbstractBlade if there is only one term
        if length(parts) == 1
            k_vectors = first(values(parts))
            if length(k_vectors) == 1
                return first(k_vectors)
            end
        end

        # Return new Multivector
        new(dim_, parts, norm_)
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

# --- Method definitions for AbstractMultivector interface functions

dim(M::Multivector) = M.dim

grades(M::Multivector; collect=true) =
    collect ?
        Base.collect(keys(M.parts)) :
        keys(M.parts)

blades(M::Multivector) = reduce(vcat, values(M.parts))

norm(M::Multivector) = M.norm

Base.getindex(M::Multivector, k::Int) =
    k in grades(M) ?  M.parts[k] : Vector{AbstractBlade}()

# --- Utility methods

convert(::Type{T}, M::Multivector) where {T<:AbstractFloat} =
    T == typeof(norm(M)) ? M : Multivector{T}(M)
