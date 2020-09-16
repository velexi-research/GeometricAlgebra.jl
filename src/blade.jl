"""
The blade.jl submodule defines Blade types.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

import Base.sign
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

For the AbstractBlade type, the norm and orientation are encoded by the `value`
of the blade. For Blades, the norm of the blade is equal to `abs(value)` and
the orientation of the blade relative to its `basis` is equal `sign(value)`.
For Scalars, the `basis` of the scalar is `1`, so the `value` of the scalar
has its expected meaning.

Methods
-------
    dim(B::AbstractBlade{T})::Int
    grade(B::AbstractBlade{T})::Int
    basis(B::AbstractBlade{T})::Union{Matrix{T}, Real}
    value(B::AbstractBlade{T})::T
    norm(B::AbstractBlade{T})::T
    sign(B::AbstractBlade{T})::Int8

Unary Operations
----------------
    -(B::AbstractBlade{T})::AbstractBlade{T}

Binary Operations
------------------
    +(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    -(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    *(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    /(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    ∧(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    ⋅(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
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
    # * `basis`: an orthonormal for the space spanned by the blade. Note that
    #   the order of the columns in `basis` defines the orientation for the
    #   unit blade represented by `basis`.
    #
    # * `value`: the signed-norm (hypervolume) of the blade. The sign
    #    of `value` indicates the orientation of the blade relative to the
    #    unit blade represented by `basis`. It is positive when the blade has
    #    the same orientation as `basis` and negative when the blade has the
    #    opposite orientation.
    dim::Int
    grade::Int
    basis::Matrix{T}
    value::T

    """
        Blade{T}(vectors::Matrix{T};
                 value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a collection of vectors stored as the columns of a
    matrix. Zero{T}() is returned when the norm of the blade is less than
    `atol`.

    By default, `vectors` determines the value of the blade. However, if
    `value` is specified, `vectors` is only used to define the subspace
    represented by the blade.

    When `value` is positive, the orientation of the blade is the same as the
    orientation of the outer product of the columns of `vectors` (taken in
    order). When `value` is negative, the orientation of the blade is the
    opposite of the orientation implied by the `vectors`.
    """
    function Blade{T}(vectors::Matrix{T};
                      value::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `value` is effectively zero
        if value != nothing && abs(value) < atol
            return Zero{T}()
        end

        # number of vectors > dimension of column space
        dims = size(vectors)
        if dims[2] > dims[1]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return Zero{T}()
            end
        end

        # --- Construct Blade

        # Preparations
        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))

        # Orthonormal basis for subspace
        basis::Matrix{T} = F.Q

        # Compute value
        value = (value == nothing) ? signed_norm : value * sign(signed_norm)

        # Return Zero if `value` is effectively zero
        if abs(value) < atol
            return Zero{T}()
        end

        # Return new Blade
        new(dims[1], dims[2], basis, value)
    end

    """
        Blade{T}(vector::Vector{T};
                 value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a single vector. Zero{T}() is returned when the
    norm of the blade is less than `atol`.

    By default, `vector` determines the value of the blade. However, if
    `value` is specified, `vector` is only used to define the subspace
    represented by the blade.

    When `value` is positive, the orientation of the blade is the same as the
    direction of `vector`. When `value` is negative, the orientation of the
    blade is the opposite of the direction of `vector`.
    """
    function Blade{T}(vector::Vector{T};
                      value::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `value` is effectively zero
        if value != nothing && abs(value) < atol
            return Zero{T}()
        end

        # --- Construct Blade

        # Compute value
        value = (value == nothing) ? LinearAlgebra.norm(vector) : value

        # Return Zero if `value` is effectively zero
        if abs(value) < atol
            return Zero{T}()
        end

        # Compute basis
        basis::Matrix{T} = reshape(vector, length(vector), 1) / abs(value)

        # Return new Blade
        new(length(vector), 1, basis, value)
    end

    """
        Blade{T}(B::Blade{T};
                 value::Union{Real, Nothing}=value(B),
                 sign::Real=sign(B), copy_basis::Bool=false)
                 where {T<:AbstractFloat}

    Construct a Blade representing the same space as `B` having a specified
    `value`. When `copy_basis` is true, the `basis` of the new Blade is a copy
    of the `basis` of the original Blade; otherwise, the `basis` of the new
    Blade is reference to the `basis` of the original Blade.
    """
    function Blade{T}(B::Blade{T};
                      value::Union{Real, Nothing}=value(B),
                      copy_basis::Bool=false) where {T<:AbstractFloat}
        # Return new Blade
        copy_basis ? new(dim(B), grade(B), copy(basis(B)), value) :
                     new(dim(B), grade(B), basis(B), value)
     end
end

"""
    Blade(vectors::Array{T};
          value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

    Blade{T}(vectors::Array{<:AbstractFloat};
             value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer};
          value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(Float64))

    Blade{T}(vectors::Array{<:Integer};
             value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

Construct a Blade from a collection of vectors represented as (1) the columns
of a matrix or (2) a single vector. Zero{T}() is returned when the norm of the
blade is less than `atol`.

By default, `vectors` determines the value (i.e., norm and orientation) of the
blade. However, if `value` is specified, `vectors` is only used to define the
subspace represented by the blade.

When `value` is positive, the orientation of the blade is the same as the
orientation of the outer product of the columns of `vectors` (taken in order).
When `value` is negative, the orientation of the blade is the opposite of the
orientation implied by the `vectors`.

When the precision is not specified, the following rules are applied to set
the precision of the Blade.

* If `vectors` is an Array of floating-point values, the precision of the
  constructed Blade is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the constructed Blade
  defaults to `Float64`.
"""
Blade(vectors::Array{T};
      value::Union{Real, Nothing}=nothing,
      atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade{T}(vectors, value=value, atol=atol)

Blade{T}(vectors::Array{<:AbstractFloat};
         value::Union{Real, Nothing}=nothing,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), value=value, atol=atol)

Blade(vectors::Array{<:Integer};
      value::Union{Real, Nothing}=nothing, atol::Real=blade_atol(Float64)) =
    Blade(convert(Array{Float64}, vectors), value=value, atol=atol)

Blade{T}(vectors::Array{<:Integer};
         value::Union{Real, Nothing}=nothing,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), value=value, atol=atol)

"""
    Blade(B::Blade{T}; value::Real=value(B), copy_basis=false)
        where {T<:AbstractFloat}

Construct a Blade representing the same space as `B` having a specified norm
and orientation relative to `B`. When `copy_basis` is true, the `basis` of the
new Blade is a copy of the `basis` of the original Blade; otherwise, the
`basis` of the new Blade is reference to the `basis` of the original Blade.
"""
Blade(B::Blade{T};
      value::Real=value(B), copy_basis=false) where {T<:AbstractFloat} =
    Blade{T}(B, value=value, copy_basis=copy_basis)

# TODO Consider adding convenience functions Blade(x::Real)


# Scalar
"""
    struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}

Scalar (0-blade) represented with the floating-point precision of type `T`.
"""
struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}
    # Fields
    # ------
    # * `value`: the value of the scalar
    value::T

    """
        Scalar{T}(value::T; atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Construct a Scalar with the specified value. Scalars with absolute value
    less than `atol` are returned as Zero{T}().
    """
    Scalar{T}(value::T; atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
        abs(value) < atol ? Zero{T}() : new(value)
end

"""
    Scalar(value::T; atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Scalar{T}(value::AbstractFloat;
              atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Scalar(value::Integer)

    Scalar{T}(value::Integer) where {T<:AbstractFloat}

Construct a Scalar with the specified value. Zero{T}() is returned when the
absolute value of `value` is (1) less than `atol` or (2) exactly equal to
zero.

When the precision is not specified, the following rules are applied to set
the precision of the Scalar.

* If `value` is a floating-point value, the precision of the constructed
  Scalar is inferred from the precision of `value`.

* If `value` is an integer, the precision of the constructed Scalar defaults
  to `Float64`.
"""
Scalar(value::T; atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Scalar{T}(value, atol=atol)

Scalar{T}(value::AbstractFloat;
          atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Scalar(convert(T, value), atol=atol)

Scalar(value::Integer) = (abs(value) == 0) ?
    Zero{Float64}() : Scalar{Float64}(convert(Float64, value))

Scalar{T}(value::Integer) where {T<:AbstractFloat} =
    (abs(value) == 0) ? Zero{T}() : Scalar{T}(convert(T, value))


# Zero
"""
    struct Zero{T<:AbstractFloat} <: AbstractScalar{T}

The additive identity 0.
"""
struct Zero{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    Zero()
    Zero(B::AbstractBlade{T}) where {T<:AbstractFloat}
    Zero(::Type{T}) where {T<:AbstractFloat}
    Zero(::Type{<:AbstractBlade})
    Zero(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat}

Return the additive identity 0. When the precision is not specified, it
defaults to `Float64`.
"""
Zero() = Zero{Float64}()
Zero(B::AbstractBlade{T}) where {T<:AbstractFloat} = Zero{T}()
Zero(::Type{T}) where {T<:AbstractFloat} = Zero{T}()
Zero(::Type{<:AbstractBlade}) = Zero{Float64}()
Zero(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat} = Zero{T}()


# One
"""
    struct One{T<:AbstractFloat} <: AbstractScalar{T}

The multiplicative identity 1.
"""
struct One{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    One()
    One(B::AbstractBlade{T}) where {T<:AbstractFloat}
    One(::Type{T}) where {T<:AbstractFloat}
    One(::Type{<:AbstractBlade})
    One(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1. When the precision is not specified, it
defaults to `Float64`.
"""
One() = One{Float64}()
One(B::AbstractBlade{T}) where {T<:AbstractFloat} = One{T}()
One(::Type{T}) where {T<:AbstractFloat} = One{T}()
One(::Type{<:AbstractBlade}) = One{Float64}()
One(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat} = One{T}()


# --- Basic Blade functions

# Exports
export dim, grade, basis, value, norm, sign

"""
    dim(B::AbstractBlade{<:AbstractFloat})

Return dimension of space that Blade blade is embedded in.
"""
dim(B::Blade) = B.dim
dim(B::AbstractScalar) = 0

"""
    grade(B::AbstractBlade{<:AbstractFloat})

Return the grade of the dimension of the space spanned by the blade.
"""
grade(B::Blade) = B.grade
grade(B::AbstractScalar) = 0

"""
    basis(B::AbstractBlade{<:AbstractFloat})

When `B` is a Blade, return an orthonormal basis for the space spanned by the
blade. When `B` is a scalar (subtype of AbstractScalar), return 1.
"""
basis(B::Blade) = B.basis
basis(B::AbstractScalar{T}) where {T<:AbstractFloat} = T(1)

"""
    value(B::AbstractBlade)::Real

Return the value of a blade. For Blades, `value(B)` is the signed norm of
the blade relative to its unit basis. For Scalars, `value(B)` is the value
of the scalar (note that the basis for Scalars is 1).
"""
value(B::Blade) = B.value
value(B::Scalar) = B.value
value(B::Zero{T}) where {T<:AbstractFloat} = T(0)
value(B::One{T}) where {T<:AbstractFloat} = T(1)

"""
    norm(B::AbstractBlade{<:AbstractFloat})

Return the norm of the blade.
"""
norm(B::Blade) = abs(value(B))
norm(B::Scalar) = abs(value(B))
norm(B::Zero) = value(B)
norm(B::One) = value(B)

"""
    sign(B::AbstractBlade)::Int8

Return the sign of a blade relative to its unit basis.
"""
sign(B::Blade)::Int8 = sign(B.value)
sign(B::Scalar)::Int8 = sign(B.value)
sign(B::Zero)::Int8 = 0
sign(B::One)::Int8 = 1


# --- Utility functions

# Exports
export blade_atol

# TODO: review numerical error in factorizations to see if a different
#       tolerance would be better.
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
