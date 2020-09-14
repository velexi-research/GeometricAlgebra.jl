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
    # * `basis`: an orthonormal for the space spanned by the blade. Note that
    #   the order of the columns in `basis` defines the orientation for the
    #   unit blade represented by `basis`.
    #
    # * `norm`: the norm (hypervolume) of the blade. It is pre-computed and
    #   cached for efficiency.
    #
    # * `sign`: the orientation of the blade relative to the unit blade
    #   represented by `basis`. It is equal to +1 when the blade has the
    #   same orientation as `basis` and is equal to -1 when the blade has
    #   the opposite orientation.
    #
    dim::Int
    grade::Int
    basis::Matrix{T}
    norm::T
    sign::Int8

    """
        Blade{T}(vectors::Matrix{T};
                 norm::Union{Real, Nothing}=nothing, sign::Real=1,
                 atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a collection of vectors stored as the columns of a
    matrix. Zero{T}() is returned when the norm of the blade is less than
    `atol`.

    By default, `vectors` determines the norm and orientation of the blade.
    However, if `norm` or `sign` are provided, `vectors` is only used to
    define the subspace represented by the blade.

    When `sign` is +1, the orientation of the blade is the same as the
    orientation of the outer product of the columns of `vectors` (taken in
    order). When `sign` is -1, the orientation of the blade is the opposite
    of the orientation implied by the `vectors`.

    When `norm` is negative, the orientation of the blade is opposite of the
    orientation determined by `sign`.
    """
    function Blade{T}(vectors::Matrix{T};
                      norm::Union{Real, Nothing}=nothing, sign::Real=1,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `norm` is effectively zero
        if norm != nothing && abs(norm) < atol
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

        # Compute norm
        if norm == nothing
            norm = abs(signed_norm)

            # blade is effectively zero
            if norm < atol
                return Zero{T}()
            end
        end

        # Compute orientation
        sign *= Base.sign(signed_norm)
        if norm < 0
            sign = -sign
        end

        # Return new Blade
        new(dims[1], dims[2], basis, norm, sign)
    end

    """
        Blade{T}(vector::Vector{T};
                 norm::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a single vector. Zero{T}() is returned when the
    norm of the blade is less than `atol`.

    By default, `vector` determines the norm and orientation of the blade.
    However, if `norm` is provided, `vector` is only used to define the
    subspace represented by the blade.

    When `sign` is +1, the orientation of the blade is the same as the
    direction of `vector`. When `sign` is -1, the orientation of the blade is
    the opposite of the direction of `vector`.

    When `norm` is negative, the orientation of the blade is opposite of the
    orientation determined `vector`.
    """
    function Blade{T}(vector::Vector{T};
                      norm::Union{Real, Nothing}=nothing, sign::Real=1,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `norm` is effectively zero
        if norm != nothing && abs(norm) < atol
            return Zero{T}()
        end

        # --- Construct Blade

        # Compute norm and orientation
        sign::Int8 = 1
        if norm == nothing
            norm = LinearAlgebra.norm(vector)

            # blade is effectively zero
            if norm < atol
                return Zero{T}()
            end
        else
            # Enforce convention that norm > 0
            if norm < 0
                sign = -1
                norm = abs(norm)
            end
        end

        # Compute basis
        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm

        # Return new Blade
        new(length(vector), 1, basis, norm, sign)
    end

    """
        Blade{T}(B::Blade{T};
                 norm::Union{Real, Nothing}=B.norm,
                 sign::Real=B.sign, copy_basis::Bool=false)
                 where {T<:AbstractFloat}

    Construct a Blade representing the same space as `B` having a specified
    norm and orientation relative to `B`. When `copy_basis` is true, the
    `basis` of the new Blade is a copy of the `basis` of the original Blade;
    otherwise, the `basis` of the new Blade is reference to the `basis` of
    the original Blade.

    When `norm` is negative, the orientation of the blade is opposite of the
    orientation determined by `sign`.
    """
    function Blade{T}(B::Blade{T};
             norm::Union{Real, Nothing}=B.norm, sign::Real=B.sign,
             copy_basis::Bool=false) where {T<:AbstractFloat}
        # `norm` < 0, so reverse sign
        if norm < 0
            norm = abs(norm)
            sign = - sign
        end

        # Return new Blade
        copy_basis ? new(dim(B), grade(B), copy(basis(B)), norm, sign) :
                     new(dim(B), grade(B), basis(B), norm, sign)
     end
end

"""
    Blade(vectors::Array{T};
          norm::Union{Real, Nothing}=nothing, sign::Real=1,
          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Blade{T}(vectors::Array{<:AbstractFloat};
             norm::Union{Real, Nothing}=nothing, sign::Real=1,
             atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer};
          norm::Union{Real, Nothing}=nothing, sign::Real=1,
          atol::Real=blade_atol(Float64))

    Blade{T}(vectors::Array{<:Integer};
             norm::Union{Real, Nothing}=nothing, sign::Real=1,
             atol::Real=blade_atol(T)) where {T<:AbstractFloat}

Construct a Blade from a collection of vectors represented as (1) the columns
of a matrix or (2) a single vector. Zero{T}() is returned when the norm of the
blade is less than `atol`.

By default, `vectors` determines the norm and orientation of the blade.
However, if `norm` or `sign` are provided, `vectors` is only used to define
the subspace represented by the blade.

When `sign` is +1, the orientation of the blade is the same as the orientation
of the outer product of the columns of `vectors` (taken in order). When `sign`
is -1, the orientation of the blade is the opposite of the orientation implied
by the `vectors`.

When `norm` is negative, the orientation of the blade is opposite of the
orientation determined by `sign`.

When the precision is not specified, the following rules are applied to set
the precision of the Blade.

* If `vectors` is an Array of floating-point values, the precision of the
  constructed Blade is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the constructed Blade
  defaults to `Float64`.
"""
Blade(vectors::Array{T};
      norm::Union{Real, Nothing}=nothing, sign::Real=1,
      atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade{T}(vectors, norm=norm, sign=sign, atol=atol)

Blade{T}(vectors::Array{<:AbstractFloat};
         norm::Union{Real, Nothing}=nothing, sign::Real=1,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), norm=norm, sign=sign, atol=atol)

Blade(vectors::Array{<:Integer};
      norm::Union{Real, Nothing}=nothing, sign::Real=1,
      atol::Real=blade_atol(Float64)) =
    Blade(convert(Array{Float64}, vectors), norm=norm, sign=sign, atol=atol)

Blade{T}(vectors::Array{<:Integer};
         norm::Union{Real, Nothing}=nothing, sign::Real=1,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), norm=norm, sign=sign, atol=atol)

"""
    Blade(B::Blade{T}; norm=B.norm, sign=B.sign, copy_basis=false)
        where {T<:AbstractFloat}

Construct a Blade representing the same space as `B` having a specified norm
and orientation relative to `B`. When `copy_basis` is true, the `basis` of the
new Blade is a copy of the `basis` of the original Blade; otherwise, the
`basis` of the new Blade is reference to the `basis` of the original Blade.
"""
Blade(B::Blade{T};
      norm=B.norm, sign=B.sign, copy_basis=false) where {T<:AbstractFloat} =
    Blade{T}(B, norm=norm, sign=sign, copy_basis=copy_basis)

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
export dim, grade, basis, norm

"""
    dim(B::AbstractBlade{<:AbstractFloat})

Return dimension of space that Blade blade is embedded in
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

Return an orthonormal for the space spanned by the blade.
"""
basis(B::Blade) = B.basis
basis(B::AbstractScalar) = nothing

"""
    norm(B::AbstractBlade{<:AbstractFloat})

Return the norm of the blade.
"""
norm(B::Blade) = B.norm
norm(B::Scalar) = abs(B.value)
norm(B::Zero) = 0
norm(B::One) = 1


# --- Basic Scalar functions

# Exports
export value

"""
    value(B::AbstractScalar{<:AbstractFloat})

Return the value of a scalar.
"""
value(B::Scalar) = B.value
value(B::Zero{T}) where {T<:AbstractFloat} = T(0)
value(B::One{T}) where {T<:AbstractFloat} = T(1)


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
