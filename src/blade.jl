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
export AbstractBlade, AbstractScalar, Blade, Scalar
export Zero, One


# AbstractBlade
"""
    abstract type AbstractBlade{T<:AbstractFloat}

Supertype for all blade types. Blades are represented with the floating-point
precision of type `T`.
"""
abstract type AbstractBlade{T<:AbstractFloat} end


# AbstractScalar
"""
    abstract type AbstractScalar{T<:AbstractFloat} <: AbstractBlade{T}

Supertype for all scalar types (i.e., 0-blades). Scalars are represented with
the floating-point precision of type `T`.
"""
abstract type AbstractScalar{T<:AbstractFloat} <: AbstractBlade{T} end


# Blade
"""
    struct Blade{T<:AbstractFloat} <: AbstractBlade{T}

Blade (having nonzero grade) represented with the floating-point precision of
type `T`.
"""
struct Blade{T<:AbstractFloat} <: AbstractBlade{T}
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
        Blade{T}(vectors::Matrix{T};
                 atol::Real=eps(T)) where {T<:AbstractFloat}

    Construct a Blade from a collection of vectors stored as the columns
    of a matrix. Blades with norm less than `atol` are returned as Zero.
    """
    function Blade{T}(vectors::Matrix{T};
                      atol::Real=eps(T)) where {T<:AbstractFloat}

        dims = size(vectors)
        if dims[1] < dims[2]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return Zero
            end
        else
            F = LinearAlgebra.qr(vectors)
            basis::Matrix{T} = F.Q
            norm::T = abs(prod(LinearAlgebra.diag(F.R)))

            if norm < atol
                return Zero
            end

            new(dims[1], dims[2], basis, norm)
        end
    end

    """
        Blade{T}(vector::Vector{T};
                 atol::Real=eps(T)) where {T<:AbstractFloat}

    Construct a Blade from a single vector. Vectors with norm less than `atol`
    are returned as Zero.
    """
    function Blade{T}(vector::Vector{T};
                      atol::Real=eps(T)) where {T<:AbstractFloat}

        norm::T = LinearAlgebra.norm(vector)

        if norm < atol
            return Zero
        end

        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm
        new(length(vector), 1, basis, norm)
    end
end

"""
    Blade(vectors::Array{T}; atol::Real=eps(T)) where {T<:AbstractFloat}

Construct a Blade from a collection of vectors represented as (1) the columns
of a matrix or (2) a single vector. Blades with norm less than `atol` are
returned as Zero.

The precision of the Blade is inferred precision from the precision of the
elements of `vectors`.
"""
Blade(vectors::Array{T}; atol::Real=eps(T)) where {T<:AbstractFloat} =
    Blade{T}(vectors, atol=atol=atol)

"""
    Blade(vectors::Array{<:Integer}; atol::Real=eps(Float64))

    Blade{T}(vectors::Array{<:Integer};
             atol::Real=eps(T)) where {T<:AbstractFloat}

Construct a Blade from a collection of vectors with integer components
represented as (1) the columns of a matrix or (2) a single vector. Blades with
norm less than `atol` are returned as Zero.

When the precision of the Blade is not explicitly specified, it defaults to
Float64.
"""
Blade(vectors::Array{<:Integer}; atol::Real=eps(Float64)) =
    Blade(convert(Array{Float64}, vectors), atol=atol)

Blade{T}(vectors::Array{<:Integer};
         atol::Real=eps(T)) where {T<:AbstractFloat} =
    Blade(vectors, atol=atol)


# Scalar
"""
    struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}

Scalar (0-blade) represented with the floating-point precision of type `T`.
"""
struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}
    value::T

    """
        Scalar{T}(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}

    Construct a Scalar with the specified value. Scalars with absolute value
    less than `atol` are returned as Zero.
    """
    function Scalar{T}(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}

        if abs(value) < atol
            return Zero
        end

        new(value)
    end
end

"""
    Scalar(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}

Construct a Scalar with the specified value. Scalars with absolute value less
than `atol` are returned as Zero.

The precision of the Blade is inferred precision from the precision of `value`.
"""
Scalar(value::T; atol::Real=eps(T)) where {T<:AbstractFloat} =
    Scalar{T}(value, atol=atol)

"""
    Scalar(value::Integer)

    Scalar{T}(value::Integer) where {T<:AbstractFloat}

Construct a Scalar with the specified value that is an Integer. When `value`
is equal to zero, Zero is return.

When the precision of the Scalar is not explicitly specified, it defaults to
Float64.
"""
Scalar(value::Integer) =
    (abs(value) == 0) ? Zero : Scalar{Float64}(convert(Float64, value))

Scalar{T}(value::Integer) where {T<:AbstractFloat} =
    (abs(value) == 0) ? Zero : Scalar{T}(convert(T, value))


# Zero
"""
    struct Zero{T<:AbstractFloat} <: AbstractScalar{T}

The additive identity 0.
"""
struct Zero{T<:AbstractFloat} <: AbstractScalar{T} end
Zero() = Zero{Float64}()

# One
"""
    struct One{T<:AbstractFloat} <: AbstractScalar{T}

The multiplicative identity 1.
"""
struct One{T<:AbstractFloat} <: AbstractScalar{T} end
One() = One{Float64}()


# --- Functions

# Exports
export dim, grade, norm, basis, inverse


# dim()
dim(B::Blade{T}) where {T<:AbstractFloat} = B.dim
dim(B::AbstractScalar{T}) where {T<:AbstractFloat} = 0


# grade()
grade(B::Blade{T}) where {T<:AbstractFloat} = B.grade
grade(B::AbstractScalar{T}) where {T<:AbstractFloat} = 0


# norm()
norm(B::Blade{T}) where {T<:AbstractFloat} = B.norm
norm(B::Scalar{T}) where {T<:AbstractFloat} = abs(B.value)
norm(B::Zero{T}) where {T<:AbstractFloat} = 0
norm(B::One{T}) where {T<:AbstractFloat} = 1


# basis()
basis(B::Blade{T}) where {T<:AbstractFloat} = B.basis
basis(B::AbstractScalar{T}) where {T<:AbstractFloat} = nothing


# inverse()
# inverse(B::Blade{T}) = Blade  # TODO
inverse(B::Scalar{T}) where {T<:AbstractFloat} = Scalar{T}(1 / B.value)
inverse(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
inverse(B::One{T}) where {T<:AbstractFloat} = One
