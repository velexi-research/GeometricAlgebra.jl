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

# Standard library
import LinearAlgebra


# --- Types

# Exports
export AbstractBlade, Scalar, Blade, Pseudoscalar

# AbstractBlade
"""
    abstract type AbstractBlade{T<:AbstractFloat}

Supertype for all blade types.

For the AbstractBlade type, the norm and orientation are encoded by the `volume`
of the blade.

Methods
-------
    dim(B)::Int
    grade(B)::Int
    basis(B)::Matrix{AbstractFloat}
    volume(B)::AbstractFloat
    norm(B)::AbstractFloat
    sign(B)::Int8

Unary Operations
----------------
    -(B)::AbstractBlade
    dual(B)::AbstractBlade
    reciprocal(B)::AbstractBlade
    reverse(B)::AbstractBlade

Binary Operations
------------------
    ∧(B, C)::AbstractBlade
    outer(B, C)::AbstractBlade

    ⋅(B, C)::AbstractBlade
    inner(B, C)::AbstractBlade

    +(B, C)::AbstractMultivector
    -(B, C)::AbstractMultivector
    *(B, C)::Union{AbstactBlade, AbstractMultivector}
    /(B, C)::AbstractMultivector

    dual(A, B)::AbstractBlade
    project(A, B)::AbstractBlade
"""
abstract type AbstractBlade{T<:AbstractFloat} end


# Scalar
"""
    struct Scalar{T<:AbstractFloat} <: AbstractBlade{T}

Scalar (0-blade) represented with the floating-point precision of type `T`. The
`basis` and `volume` of a Scalar are `1` and the value of the Scalar,
respectively.
"""
struct Scalar{T<:AbstractFloat} <: AbstractBlade{T}
    #=
      Fields
      ------
      * `value`: the value of the scalar
    =#
    value::T

    """
        Scalar{T}(value::T) where {T<:AbstractFloat}

    Construct a Scalar with the specified value.
    """
    Scalar{T}(value::T) where {T<:AbstractFloat} = new(value)
end

"""
    Scalar(value::T) where {T<:AbstractFloat}
    Scalar{T}(value::AbstractFloat) where {T<:AbstractFloat}
    Scalar(value::Integer)
    Scalar{T}(value::Integer) where {T<:AbstractFloat}

Construct a Scalar with the specified value.

When the precision is not specified, the following rules are applied to set
the precision of the Scalar.

* If `value` is a floating-point value, the precision of the Scalar is inferred
  from the precision of `value`.

* If `value` is an integer, the precision of the Scalar defaults to `Float64`.
"""
Scalar(value::T) where {T<:AbstractFloat} = Scalar{T}(value)
Scalar{T}(value::AbstractFloat) where {T<:AbstractFloat} = Scalar(T(value))
Scalar(value::Integer) = Scalar(Float64(value))
Scalar{T}(value::Integer) where {T<:AbstractFloat} = Scalar(T(value))


# Blade
"""
    struct Blade{T<:AbstractFloat} <: AbstractBlade{T}

Blade (having nonzero grade) represented with the floating-point precision of
type `T`. The norm and orientation of a Blade are encoded by its `volume`. The
norm of a Blade is equal to `abs(volume)` and the orientation of a Blade
relative to its `basis` is equal to `sign(volume)`.
"""
struct Blade{T<:AbstractFloat} <: AbstractBlade{T}
    #=
     Fields
      ------
      * `dim`: the dimension of the space that the blade is embedded in

      * `grade`: the dimension of the space spanned by the blade

      * `basis`: an orthonormal for the space spanned by the blade. Note that
        the order of the columns in `basis` defines the orientation for the
        unit blade represented by `basis`.

      * `volume`: the signed-norm (hypervolume) of the blade. The sign
         of `volume` indicates the orientation of the blade relative to the
         unit blade represented by `basis`. It is positive when the blade has
         the same orientation as `basis` and negative when the blade has the
         opposite orientation.
    =#
    dim::Int
    grade::Int
    basis::Matrix{T}
    volume::T

    """
        Blade{T}(dim::Int, grade::Int,
                 basis::Matrix{T},
                 volume::Real;
                 atol::Real=blade_atol(T),
                 enforce_constraints::Bool=true,
                 copy_basis::Bool=true)
            where {T<:AbstractFloat}

    Construct a Blade from the specified data field values. If the absolute
    value of `volume` is less than `atol`, a Scalar representing zero is
    returned. If the `dim` and `grade` are equal, a Pseudoscalar is returned.

    When `enforce_constraints` is true, constraints are enforced. When
    `copy_basis` is true, the basis of the new Blade is a copy of `basis`;
    otherwise, the basis of the new Blade is a reference to `basis`.

    Note: this inner constructor intended primarily for internal use by
    other inner constructors to enforce constraints.
    """
    function Blade{T}(dim::Int,
                      grade::Int,
                      basis::Matrix{T},
                      volume::Real;
                      atol::Real=blade_atol(T),
                      enforce_constraints::Bool=true,
                      copy_basis::Bool=true) where {T<:AbstractFloat}

        # --- Enforce constraints

        if enforce_constraints
            basis_size = size(basis)

            if dim != basis_size[1]
                message = string("`dim` not equal to the number of ",
                                 "rows in `basis`")
                throw(DimensionMismatch(message))
            end

            if grade != basis_size[2]
                message = string("`grade` not equal to the number of ",
                                 "columns in `basis`")
                throw(DimensionMismatch(message))
            end

            for i in basis_size[2]
                if LinearAlgebra.norm(basis[:, i]) ≉ 1
                    error("columns of `basis` not normalized")
                end
            end
        end

        # --- Handle edge cases

        # `volume` is effectively zero
        if abs(volume) < atol
            return zero(Blade{T})
        end

        # Return a Pseudoscalar if the grade of the blade is equal to the
        # dimension of the embedding space.
        if grade == dim
            return Pseudoscalar{T}(dim, volume)
        end

        # Return new Blade
        copy_basis ?
        new(dim, grade, copy(basis), volume) : new(dim, grade, basis, volume)
    end

    """
        Blade{T}(vectors::Matrix{T};
                 volume::Union{Real, Nothing}=nothing,
                 atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a collection of vectors stored as the columns of a
    matrix. If the norm the blade is less than `atol`, a Scalar representing
    zero is returned. If the grade of the blade is equal to the dimension of
    the space that the blade is embedded in, a Pseudoscalar is returned.

    By default, `vectors` determines the volume of the blade. However, if
    `volume` is specified, `vectors` is only used to define the subspace
    (including orientation) represented by the blade.

    When `volume` is positive, the orientation of the blade is the same as the
    orientation of the outer product of the columns of `vectors` (taken in
    order). When `volume` is negative, the orientation of the blade is the
    opposite of the orientation implied by the `vectors`.
    """
    function Blade{T}(vectors::Matrix{T};
                      volume::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `volume` is effectively zero
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
        end

        # number of vectors > dimension of column space
        dims = size(vectors)
        if dims[2] > dims[1]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return zero(Blade{T})
            end
        end

        # --- Construct Blade

        # Preparations
        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))

        # Return zero if norm is below tolerance
        if abs(signed_norm) < atol
            return zero(Blade{T})
        end

        # Compute orthonormal basis for subspace
        basis::Matrix{T} = F.Q

        # Compute volume
        volume = (volume == nothing) ?  signed_norm : volume * sign(signed_norm)

        # Return new Blade
        Blade{T}(dims[1], dims[2], basis, volume,
                 atol=atol, enforce_constraints=false, copy_basis=false)
    end

    """
        Blade{T}(v::Vector{T};
                 volume::Union{Real, Nothing}=nothing,
                 atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a single vector. A Scalar representing zero is
    returned when the norm of the blade is less than `atol`.

    By default, `v` determines the volume of the blade. However, if `volume`
    is specified, `v` is only used to define the subspace represented by the
    blade.

    When `volume` is positive, the orientation of the blade is the same as the
    direction of `v`. When `volume` is negative, the orientation of the blade
    is the opposite of the direction of `v`.
    """
    function Blade{T}(v::Vector{T};
                      volume::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `volume` is effectively zero
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
        end

        # --- Construct Blade

        norm_v = LinearAlgebra.norm(v)

        # Return zero if norm is below tolerance
        if abs(norm_v) < atol
            return zero(Blade{T})
        end

        # Compute basis
        basis::Matrix{T} = reshape(v, length(v), 1) / norm_v

        # Compute volume
        volume = (volume == nothing) ? norm_v : volume

        # Return new Blade
        Blade{T}(length(v), 1, basis, volume, atol=atol,
                 enforce_constraints=false, copy_basis=false)
    end

    """
        Blade{T}(B::Blade{T};
                 volume::Real=volume(B),
                 atol::Real=blade_atol(T),
                 copy_basis::Bool=false)
            where {T<:AbstractFloat}

    Construct a Blade representing the same space as `B` having a specified
    oriented `volume` relative to `B`. A Scalar representing zero is returned
    if the absolute value of `volume` is less than `atol`,

    When `copy_basis` is true, the `basis` of the new Blade is a copy
    of the `basis` of the original Blade; otherwise, the `basis` of the new
    Blade is reference to the `basis` of the original Blade.
    """
    Blade{T}(B::Blade{T};
             volume::Real=volume(B),
             atol::Real=blade_atol(T),
             copy_basis::Bool=false) where {T<:AbstractFloat} =

        Blade{T}(dim(B), grade(B), basis(B), volume,
                 atol=atol, enforce_constraints=false, copy_basis=copy_basis)

    """
        Blade{T}(B::Blade{S};
                 volume::Real=volume(B),
                 atol::Real=blade_atol(T),
                 copy_basis::Bool=false)
            where {T<:AbstractFloat, S<:AbstractFloat}

    Convert the floating-point precision of a Blade.

    If `volume` is specified, the oriented volume of the blade of the new blade
    (relative to `B`) is set to `volume`. A Scalar representing zero is
    returned if the absolute value of `volume` is less than `atol`.

    When `copy_basis` is true, the `basis` of the new Blade is a copy of the
    `basis` of the original Blade; otherwise, the `basis` of the new Blade is
    reference to the `basis` of the original Blade.
    """
    Blade{T}(B::Blade{S};
             volume::Real=volume(B),
             atol::Real=blade_atol(T),
             copy_basis::Bool=false) where {T<:AbstractFloat,
                                            S<:AbstractFloat} =
        Blade{T}(dim(B), grade(B), convert(Array{T}, basis(B)), volume,
                 atol=atol, enforce_constraints=false, copy_basis=false)
end

"""
    Blade(vectors::Array{T};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

    Blade{T}(vectors::Array{<:AbstractFloat};
             volume::Union{Real, Nothing}=nothing,
             atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(Float64))

    Blade{T}(vectors::Array{<:Integer};
             volume::Union{Real, Nothing}=nothing,
             atol::Real=blade_atol(T))
        where {T<:AbstractFloat}

Construct a Blade from a collection of vectors represented as (1) the columns
of a matrix or (2) a single vector. A Scalar representing zero is returned when
the norm of the blade is less than `atol`.

By default, `vectors` determines the volume (i.e., norm and orientation) of the
blade. However, if `volume` is specified, `vectors` is only used to define the
subspace represented by the blade.

When `volume` is positive, the orientation of the blade is the same as the
orientation of the outer product of the columns of `vectors` (taken in order).
When `volume` is negative, the orientation of the blade is the opposite of the
orientation implied by the `vectors`.

When the precision is not specified, the following rules are applied to set
the precision of the Blade.

* If `vectors` is an Array of floating-point values, the precision of the Blade
  is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the Blade defaults to
 `Float64`.
"""
Blade(vectors::Array{T};
      volume::Union{Real, Nothing}=nothing,
      atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade{T}(vectors, volume=volume, atol=atol)

Blade{T}(vectors::Array{<:AbstractFloat};
         volume::Union{Real, Nothing}=nothing,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), volume=volume, atol=atol)

Blade(vectors::Array{<:Integer};
      volume::Union{Real, Nothing}=nothing,
      atol::Real=blade_atol(Float64)) =
    Blade(convert(Array{Float64}, vectors), volume=volume, atol=atol)

Blade{T}(vectors::Array{<:Integer};
         volume::Union{Real, Nothing}=nothing,
         atol::Real=blade_atol(T)) where {T<:AbstractFloat} =
    Blade(convert(Array{T}, vectors), volume=volume, atol=atol)

"""
    Blade(B::Blade{T};
          volume::Real=volume(B),
          atol::Real=blade_atol(T),
          copy_basis=false) where {T<:AbstractFloat}

Copy constructor. Construct a Blade representing the same space as `B` having
a specified oriented volume relative to `B`. A Scalar representing zero is
returned if the absolute value of `volume` is less than `atol`.

When `copy_basis` is true, the `basis` of the new Blade is a copy of the
`basis` of the original Blade; otherwise, the `basis` of the new Blade is
reference to the `basis` of the original Blade.
"""
Blade(B::Blade{T};
      volume::Real=volume(B),
      atol::Real=blade_atol(T),
      copy_basis=false) where {T<:AbstractFloat} =
    Blade{T}(B, volume=volume, atol=atol, copy_basis=copy_basis)

"""
    Blade(x::T) where {T<:AbstractFloat}
    Blade{T}(x::AbstractFloat) where {T<:AbstractFloat}
    Blade(x::Integer)
    Blade{T}(x::Integer) where {T<:AbstractFloat}

Convenience constructors for constructing Scalars.
"""
Blade(x::T) where {T<:AbstractFloat} = Scalar(x)
Blade{T}(x::AbstractFloat) where {T<:AbstractFloat} = Scalar{T}(x)
Blade(x::Integer) = Scalar(x)
Blade{T}(x::Integer) where {T<:AbstractFloat} = Scalar{T}(x)

"""
    Blade(dim::Integer, x::T) where {T<:AbstractFloat}
    Blade{T}(dim::Integer, x::AbstractFloat) where {T<:AbstractFloat}
    Blade(dim::Integer, x::Integer)
    Blade{T}(dim::Integer, x::Integer) where {T<:AbstractFloat}

Convenience constructors for constructing Pseudoscalars.
"""
Blade(dim::Integer, x::T) where {T<:AbstractFloat} = Pseudoscalar(dim, x)
Blade{T}(dim::Integer, x::AbstractFloat) where {T<:AbstractFloat} =
    Pseudoscalar{T}(dim, x)
Blade(dim::Integer, x::Integer) = Pseudoscalar(dim, x)
Blade{T}(dim::Integer, x::Integer) where {T<:AbstractFloat} =
    Pseudoscalar{T}(dim, x)


# Pseudoscalar
"""
    struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade

Pseudoscalar ((n-1)-blade) represented with the floating-point precision of
type `T`. The `basis` for a Pseudoscalar is the standard basis for ``R^n``.
The norm and orientation of a Pseudoscalar are encoded in its `value`. The
norm of a Pseudoscalar is equal to `abs(value)` and the orientation of a
Pseudoscalar relative to the standard basis for ``R^n`` is equal to
`sign(value)`.
"""
struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade{T}
    #=
      Fields
      ------
      * `dim`: the dimension of the space that the blade is embedded in

      * `value`: the value of the pseudoscalar
    =#
    dim::Int
    value::T

    """
        Pseudoscalar{T}(value::T) where {T<:AbstractFloat}

    Construct a Pseudoscalar with for a geometric algebra in `dim` dimensions
    with the specified value.
    """
    Pseudoscalar{T}(dim::Integer, value::T) where {T<:AbstractFloat} =
        new(dim, value)
end

"""
    Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat}
    Pseudoscalar{T}(dim::Integer, value::AbstractFloat) where {T<:AbstractFloat}
    Pseudoscalar(dim::Integer, value::Integer)
    Pseudoscalar{T}(dim::Integer, value::Integer) where {T<:AbstractFloat}

Construct a Pseudoscalar with for a geometric algebra in `dim` dimensions with
the specified value.

When the precision is not specified, the following rules are applied to set
the precision of the Pseudoscalar.

* If `value` is a floating-point value, the precision of the Pseudoscalar is
  inferred from the precision of `value`.

* If `value` is an integer, the precision of the Pseudoscalar defaults to
  `Float64`.
"""
Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat} =
    Pseudoscalar{T}(dim, value)

Pseudoscalar{T}(dim::Integer, value::AbstractFloat) where {T<:AbstractFloat} =
    Pseudoscalar(dim, T(value))

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))

Pseudoscalar{T}(dim::Integer, value::Integer) where {T<:AbstractFloat} =
    Pseudoscalar(dim, T(value))

"""
    Pseudoscalar(B::Pseudoscalar{T};
                 value::Real=value(B)) where {T<:AbstractFloat}

Copy constructor. Construct a Pseudoscalar representing the same space as
`B` having the specified value.
"""
Pseudoscalar(B::Pseudoscalar{T};
             value::Real=value(B)) where {T<:AbstractFloat} =
    Pseudoscalar{T}(dim(B), value)


# --- Special number functions

# Imports
import Base.zero, Base.one

"""
    zero(B::AbstractBlade)
    zero(::Type{Scalar{T}}) where {T<:AbstractFloat}
    zero(::Type{Blade{T}}) where {T<:AbstractFloat}
    zero(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(B::AbstractBlade) = zero(typeof(B))
zero(::Type{Scalar{T}}) where {T<:AbstractFloat} = Scalar{T}(0)
zero(::Type{Blade{T}}) where {T<:AbstractFloat} = zero(Scalar{T})
zero(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat} = zero(Scalar{T})

"""
    one(B::AbstractBlade)
    one(::Type{Scalar{T}}) where {T<:AbstractFloat}
    one(::Type{Blade{T}}) where {T<:AbstractFloat}
    one(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(B::AbstractBlade) = one(typeof(B))
one(::Type{Scalar{T}}) where {T<:AbstractFloat} = Scalar{T}(1)
one(::Type{Blade{T}}) where {T<:AbstractFloat} = one(Scalar{T})
one(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat} = one(Scalar{T})


# --- Core AbstractBlade functions

# Imports
import Base.sign
import LinearAlgebra.norm

# Exports
export dim, grade, basis, volume, norm

"""
    dim(B::AbstractBlade)::Integer

Return dimension of space that `B` is embedded in.
"""
dim(B::Scalar) = 0
dim(B::Blade) = B.dim
dim(B::Pseudoscalar) = B.dim

"""
    grade(B::AbstractBlade)::Integer

Return the grade of the dimension of the space spanned by `B`.
"""
grade(B::Scalar) = 0
grade(B::Blade) = B.grade
grade(B::Pseudoscalar) = B.dim

"""
    basis(B::AbstractBlade)

When `B` is a Blade, return an orthonormal basis for the space spanned by the
blade. When `B` is a Scalar, return 1. When `B` is a Pseudoscalar, return
LinearAlgebra.I.
"""
basis(B::Scalar) = 1
basis(B::Blade) = B.basis
basis(B::Pseudoscalar) = LinearAlgebra.I

"""
    volume(B::AbstractBlade)::Real

Return the volume of `B`. For Blades, `volume(B)` is the signed norm of the
blade relative to its unit basis. For Scalars, `volume(B)` is the value of the
scalar (note that the basis for Scalars is 1).
"""
volume(B::Scalar) = B.value
volume(B::Blade) = B.volume
volume(B::Pseudoscalar) = B.value

"""
    norm(B::AbstractBlade)::Real

Return the norm of `B`.
"""
norm(B::AbstractBlade) = abs(volume(B))

"""
    sign(B::AbstractBlade)::Int8

Return the sign of `B` relative to its unit basis.
"""
sign(B::AbstractBlade)::Int8 = sign(volume(B))


# --- Core Scalar functions

# Exports
export value

"""
    value(B::Scalar)::Real
    value(B::Pseudoscalar)::Real

Return the value of `B`.
"""
value(B::Scalar) = B.value
value(B::Pseudoscalar) = B.value


# --- Utility functions

# Imports
import Base.convert

# Exports
export blade_atol

"""
    convert(::Type{S}, B::Scalar) where {T<:AbstractFloat, S<:Scalar{T}}

Convert Scalar to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Scalar) where {T<:AbstractFloat, S<:AbstractBlade{T}} =
    T == typeof(volume(B)) ? B : Scalar{T}(value(B))

"""
    convert(::Type{S}, B::Blade) where {T<:AbstractFloat, S<:Blade{T}}

Convert Blade to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Blade) where {T<:AbstractFloat, S<:AbstractBlade{T}} =
    T == typeof(volume(B)) ? B : Blade{T}(B)

"""
    convert(::Type{S}, B::Pseudoscalar)
        where {T<:AbstractFloat, S<:Pseudoscalar{T}}

Convert Pseudoscalar to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Pseudoscalar) where {T<:AbstractFloat,
                                           S<:AbstractBlade{T}} =
    T == typeof(volume(B)) ? B : Pseudoscalar{T}(dim(B), value(B))

#=
 TODO: review numerical error in factorizations to see if a different
       tolerance would be better.
=#
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
