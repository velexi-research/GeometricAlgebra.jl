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

import Base.:(==), Base.:(≈), Base.:(-)
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
    sign::Int16

    """
        Blade{T}(vectors::Matrix{T}; atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a collection of vectors stored as the columns
    of a matrix. Blades with norm less than `atol` are returned as Zero{T}().
    """
    function Blade{T}(vectors::Matrix{T};
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        dims = size(vectors)
        if dims[1] < dims[2]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return Zero{T}()
            end
        else
            F = LinearAlgebra.qr(vectors)
            basis::Matrix{T} = F.Q
            norm::T = abs(prod(LinearAlgebra.diag(F.R)))

            if norm < atol
                return Zero{T}()
            end

            new(dims[1], dims[2], basis, norm, 1)
        end
    end

    """
        Blade{T}(vector::Vector{T}; atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a single vector. Vectors with norm less than `atol`
    are returned as Zero{T}().
    """
    function Blade{T}(vector::Vector{T};
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        norm::T = LinearAlgebra.norm(vector)

        if norm < atol
            return Zero{T}()
        end

        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm
        new(length(vector), 1, basis, norm, 1)
    end

    """
        Blade{T}(B::AbstractBlade{T};
                 norm=B.norm, sign=B.sign, copy_basis=false)
                 where {T<:AbstractFloat}

    Construct a Blade representing the same space as `B` having a specified
    norm and orientation relative to `B`. When `copy_basis` is true, the
    `basis` of the new Blade is a copy of the `basis` of the original Blade;
    otherwise, the `basis` of the new Blade is reference to the `basis` of
    the original Blade.
    """
    Blade{T}(B::AbstractBlade{T};
             norm=B.norm, sign=B.sign,
             copy_basis=false) where {T<:AbstractFloat} =
        copy_basis ? new(dim(B), grade(B), copy(basis(B)), norm, sign) :
                     new(dim(B), grade(B), basis(B), norm, sign)
end

"""
    Blade(vectors::Array{T}; atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Blade{T}(vectors::Array{<:AbstractFloat}; atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer}; atol::Real=blade_atol(Float64))

    Blade{T}(vectors::Array{<:Integer}; atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

Construct a Blade from a collection of vectors represented as (1) the columns
of a matrix or (2) a single vector. Zero{T}() is returned when the norm of the
blade is less than `atol`.

When the precision is not specified, the following rules are applied to set
the precision of the Blade.

* If `vectors` is an Array of floating-point values, the precision of the
  constructed Blade is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the constructed Blade
  defaults to `Float64`.
"""
Blade(vectors::Array{T}; atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade{T}(vectors, atol=atol)

Blade{T}(vectors::Array{<:AbstractFloat};
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), atol=atol)

Blade(vectors::Array{<:Integer}; atol::Real=blade_atol(Float64)) =
    Blade(convert(Array{Float64}, vectors), atol=atol)

Blade{T}(vectors::Array{<:Integer};
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), atol=atol)

"""
    Blade(B::AbstractBlade{T}; norm=B.norm, sign=B.sign, copy_basis=false)
        where {T<:AbstractFloat}

Construct a Blade representing the same space as `B` having a specified norm
and orientation relative to `B`. When `copy_basis` is true, the `basis` of the
new Blade is a copy of the `basis` of the original Blade; otherwise, the
`basis` of the new Blade is reference to the `basis` of the original Blade.
"""
Blade(B::AbstractBlade{T};
      norm=B.norm, sign=B.sign, copy_basis=false) where {T<:AbstractFloat} =
    Blade{T}(B, norm=norm, sign=sign, copy_basis=copy_basis)


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


# --- Basic functions

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


# --- Comparison functions

"""
    ==(B1::AbstractBlade{<:AbstractFloat}, B2::AbstractBlade{<:AbstractFloat})

Return true if B1 and B2 are equal; otherwise, return false.
"""
==(B1::Blade, B2::Blade) =  # TODO: fix basis and orientation comparison
    B1.dim == B2.dim && B1.grade == B2.grade &&
    B1.norm == B2.norm &&
    B1.basis == B2.basis && B1.sign == B2.sign

==(B1::Scalar, B2::Scalar) = B1.value == B2.value
==(B::Scalar, x::Real) = (x == B.value)
==(x::Real, B::Scalar) = (B == x)

==(B::Zero, x::Real) = (x == 0)
==(x::Real, B::Zero) = (B == 0)
==(B1::Zero, B2::Zero) = true

==(B::One, x::Real) = (x == 1)
==(x::Real, B::One) = (B == 1)
==(B1::One, B2::One) = true

"""
    ≈(B1::AbstractBlade{<:AbstractFloat}, B2::AbstractBlade{<:AbstractFloat})

Return true if B1 and B2 are approximatly equal; otherwise, return false.
"""
≈(B1::Blade{T1}, B2::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    B1.sign == B2.sign &&
    ≈(B1.norm, B2.norm, atol=atol, rtol=rtol) &&
    true  # TODO: add check that basis represent the same space

≈(B1::Scalar{T1}, B2::Scalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(B1.value, B2.value, atol=atol, rtol=rtol)

≈(B::Scalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, B.value, rtol=rtol, atol=atol)
≈(x::Real, B::Scalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, rtol=rtol, atol=atol)


# --- Basic operations

# Exports
export inverse

"""
    -(B::AbstractBlade)

Return the additive inverse of `B`.
"""
-(B::Blade{<:AbstractFloat}) = Blade(B, norm=norm(B), sign=-1)
-(B::Scalar) = Scalar(-B.value)

"""
    inverse(B::AbstractBlade{<:AbstractFloat})

Return the multiplicative inverse of `B`.
"""
inverse(B::Blade{<:AbstractFloat}) =
    mod(grade(B), 4) < 2 ?
    Blade(B, norm=1 / norm(B), sign=1) : Blade(B, norm=1 / norm(B), sign=-1)

inverse(B::Scalar{<:AbstractFloat}) = Scalar(1 / B.value)
inverse(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
inverse(B::One{T}) where {T<:AbstractFloat} = One{T}()


# --- Helper functions

# Exports
export blade_atol

# TODO: review numerical error in factorizations to see if a different
#       tolerance would be better.
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
