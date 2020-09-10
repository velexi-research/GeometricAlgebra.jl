"""
The blade.jl submodule defines Blade types.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

import LinearAlgebra


# --- Types

# Exports
export Blade, ZeroBlade, NullBlade


# AbstractBlade type
abstract type AbstractBlade end


# Blade type
"""
    struct Blade{T<:AbstractFloat}

The Blade type represents a blade that is stored with the floating-point
precision of type `T`.
"""
struct Blade{T<:AbstractFloat} <: AbstractBlade
    # Fields
    # ------
    # * `dim`: the dimension of the space that the blade is embedded in
    #
    # * `grade`: the dimension of the space spanned by the blade
    #
    # * `basis`: an orthonormal for the space spanned by the blade
    #
    # * `norm`: the norm of the blade. It is pre-computed and cached for
    #   efficiency.
    #
    # Notes
    # -----
    # * The orientation of the blade is implicit in the order of the columns
    #   in the `basis` matrix.
    dim::Int
    grade::Int
    basis::Matrix{T}
    norm::T

    """
        Blade{T}(vectors::Array{T,2}, atol::AbstractFloat=eps(T))

    Construct a Blade from a collection of vectors stored as the columns
    of a 2-dimensional array.
    """
    function Blade{T}(vectors::Array{T,2},
                      atol::AbstractFloat=eps(T)) where {T<:AbstractFloat}

        dims = size(vectors)
        if dims[1] < dims[2]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return NullBlade
            end
        else
            F = LinearAlgebra.qr(vectors)
            basis::Matrix{T} = F.Q
            norm::T = abs(prod(LinearAlgebra.diag(F.R)))

            if norm < atol
                return NullBlade
            end

            new(dims[1], dims[2], basis, norm)
        end
    end

    """
        Blade{T}(vectors::Array{T,1})

    Construct a Blade from a single vector (1-dimensional Array)
    """
    function Blade{T}(vector::Array{T,1}) where {T<:AbstractFloat}
        norm::T = LinearAlgebra.norm(vector)
        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm
        new(length(vector), 1, basis, norm)
    end
end

"""
    Blade(vectors)

Construct a Blade from a collection of vectors stored as the columns of a
2-dimensional array of floating-point values.

The precision of the Blade is inferred precision from the precision of the
`vectors` Array.
"""
Blade(vectors::Array{T}) where {T<:AbstractFloat} = Blade{T}(vectors)

"""
    Blade(vectors::Array{<:Integer})
    Blade{T}(vectors::Array{<:Integer}) where {T<:AbstractFloat}

Construct a Blade from a collection of vectors stored as the columns of a
2-dimensional array of integer values.

When the precision of the Blade is not explicitly specified, it defaults to
Float64.
"""
Blade(vectors::Array{<:Integer}) = Blade(convert(Array{Float64}, vectors))
Blade{T}(vectors::Array{<:Integer}) where {T<:AbstractFloat} = Blade(vectors)


# ZeroBlade type
"""
    struct ZeroBlade{T<:AbstractFloat} <: Blade{T}

The ZeroBlade type represents a 0-blade with a value that is stored with the
floating-point precision of type `T`.
"""
struct ZeroBlade{T<:AbstractFloat} <: AbstractBlade
    value::T
end

"""
    ZeroBlade(value::Integer)
    ZeroBlade{T}(value::Integer) where {T<:AbstractFloat}

Construct a ZeroBlade from a scalar value that is an Integer type.

When the precision of the Blade is not explicitly specified, it defaults to
Float64.
"""
ZeroBlade(value::Integer) = ZeroBlade{Float64}(convert(Float64, value))
ZeroBlade{T}(value::Integer) where {T<:AbstractFloat} =
    ZeroBlade{T}(convert(T, value))


# NullBlade type
"""
    struct NullBlade{T<:AbstractFloat} <: Blade{T}

The NullBlade type represents 0 - the result when a blade is constructed from
a collection of linearly dependent vectors.
"""
struct NullBlade <: AbstractBlade end


# --- Functions

# Exports
export norm

# dim()
dim(B::Blade{T}) where {T<:AbstractFloat} = B.dim

# grade()
grade(B::Blade{T}) where {T<:AbstractFloat} = B.grade
grade(B::ZeroBlade{T}) where {T<:AbstractFloat} = 1

# norm()
norm(B::Blade{T}) where {T<:AbstractFloat} = B.norm
norm(B::ZeroBlade{T}) where {T<:AbstractFloat} = B.value
norm(B::NullBlade) = 0


# basis()
basis(B::Blade{T}) where {T<:AbstractFloat} = B.basis
dim(B::ZeroBlade{T}) where {T<:AbstractFloat} = nothing
basis(B::NullBlade) = nothing


# inverse()
# inverse(B::Blade{T}) = Blade  # TODO
inverse(B::ZeroBlade{T}) where {T<:AbstractFloat} =
    B.value != 0 ? 1 / B.value : NaN

inverse(B::NullBlade) = NaN
