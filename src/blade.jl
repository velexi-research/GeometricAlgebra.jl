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
import Base.convert
import LinearAlgebra


# --- Types

# Exports
export AbstractBlade, Blade, Scalar, Pseudoscalar

# AbstractBlade
"""
    abstract type AbstractBlade

Supertype for all blade types.

For the AbstractBlade type, the norm and orientation are encoded by the `volume`
of the blade. For Blades, the norm of the blade is equal to `abs(volume)` and
the orientation of the blade relative to its `basis` is equal `sign(volume)`.
For Scalars, the `basis` and `volume` of the scalar are `1` and the value
of the scalar, respectively.

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
abstract type AbstractBlade end


# Blade
"""
    struct Blade{T<:AbstractFloat} <: AbstractBlade

Blade (having nonzero grade) represented with the floating-point precision of
type `T`.
"""
struct Blade{T<:AbstractFloat} <: AbstractBlade
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

    Construct a Blade from the specified data field values. A Scalar
    representing zero is returned when the absolute value of `volume` is less
    than `atol`.

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
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
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
    matrix. A Scalar representing zero is returned when the absolute value of
    `volume` is less than `atol`.

    By default, `vectors` determines the volume of the blade. However, if
    `volume` is specified, `vectors` is only used to define the subspace
    represented by the blade.

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

        # Orthonormal basis for subspace
        basis::Matrix{T} = F.Q

        # Compute volume
        volume = (volume == nothing) ? signed_norm : volume * sign(signed_norm)

        # Return new Blade
        Blade{T}(dims[1], dims[2], basis, volume,
                 atol=atol, enforce_constraints=false, copy_basis=false)
    end

    """
        Blade{T}(vector::Vector{T};
                 volume::Union{Real, Nothing}=nothing,
                 atol::Real=blade_atol(T))
            where {T<:AbstractFloat}

    Construct a Blade from a single vector. A Scalar representing zero is
    returned when the norm of the blade is less than `atol`.

    By default, `vector` determines the volume of the blade. However, if
    `volume` is specified, `vector` is only used to define the subspace
    represented by the blade.

    When `volume` is positive, the orientation of the blade is the same as the
    direction of `vector`. When `volume` is negative, the orientation of the
    blade is the opposite of the direction of `vector`.
    """
    function Blade{T}(vector::Vector{T};
                      volume::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `volume` is effectively zero
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
        end

        # --- Construct Blade

        # Compute basis
        norm_vector = LinearAlgebra.norm(vector)
        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm_vector

        # Set volume
        volume = (volume == nothing) ? norm_vector : volume

        # Return new Blade
        Blade{T}(length(vector), 1, basis, volume, atol=atol,
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
    function Blade{T}(B::Blade{S};
                      volume::Real=volume(B),
                      atol::Real=blade_atol(T),
                      copy_basis::Bool=false) where {T<:AbstractFloat,
                                                     S<:AbstractFloat}
        # Handle special cases
        if T == S
            return B{T}(B, volume=volume, atol=atol,
                        enforce_constraints=false, copy_basis=copy_basis)
        end

        # Return new Blade
        Blade{T}(dim(B), grade(B), convert(Array{T}, basis(B)), volume,
                 atol=atol, enforce_constraints=false, copy_basis=false)
    end
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

* If `vectors` is an Array of floating-point values, the precision of the
  constructed Blade is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the constructed Blade
  defaults to `Float64`.
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


# Scalar
"""
    struct Scalar{T<:AbstractFloat} <: AbstractBlade

Scalar (0-blade) represented with the floating-point precision of type `T`.
"""
struct Scalar{T<:AbstractFloat} <: AbstractBlade
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

* If `value` is a floating-point value, the precision of the constructed
  Scalar is inferred from the precision of `value`.

* If `value` is an integer, the precision of the constructed Scalar defaults
  to `Float64`.
"""
Scalar(value::T) where {T<:AbstractFloat} = Scalar{T}(value)
Scalar{T}(value::AbstractFloat) where {T<:AbstractFloat} = Scalar(T(value))
Scalar(value::Integer) = Scalar(Float64(value))
Scalar{T}(value::Integer) where {T<:AbstractFloat} = Scalar(T(value))


# Pseudoscalar
"""
    struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade

Pseudoscalar ((n-1)-blade) represented with the floating-point precision of
type `T`.
"""
struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade
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

* If `value` is a floating-point value, the precision of the constructed
  Pseudoscalar is inferred from the precision of `value`.

* If `value` is an integer, the precision of the constructed Pseudoscalar
  defaults to `Float64`.
"""
Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat} =
    Pseudoscalar{T}(dim, value)

Pseudoscalar{T}(dim::Integer, value::AbstractFloat) where {T<:AbstractFloat} =
    Pseudoscalar(dim, T(value))

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))

Pseudoscalar{T}(dim::Integer, value::Integer) where {T<:AbstractFloat} =
    Pseudoscalar(dim, T(value))


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

# Exports
export dim, grade, basis, volume, norm, sign

"""
    dim(B::AbstractBlade)::Integer

Return dimension of space that Blade blade is embedded in.
"""
dim(B::Blade) = B.dim
dim(B::Scalar) = 0
dim(B::Pseudoscalar) = B.dim

"""
    grade(B::AbstractBlade)::Integer

Return the grade of the dimension of the space spanned by the blade.
"""
grade(B::Blade) = B.grade
grade(B::Scalar) = 0
grade(B::Pseudoscalar) = B.dim

"""
    basis(B::AbstractBlade)

When `B` is a Blade, return an orthonormal basis for the space spanned by the
blade. When `B` is a Scalar, return 1.
"""
basis(B::Blade) = B.basis
basis(B::Scalar) = 1
basis(B::Pseudoscalar) = LinearAlgebra.I

"""
    volume(B::AbstractBlade)::Real

Return the volume of a blade. For Blades, `volume(B)` is the signed norm of
the blade relative to its unit basis. For Scalars, `volume(B)` is the value
of the scalar (note that the basis for Scalars is 1).
"""
volume(B::Blade) = B.volume
volume(B::Scalar) = B.value
volume(B::Pseudoscalar) = B.value

"""
    norm(B::AbstractBlade)::Real

Return the norm of the blade.
"""
norm(B::Blade) = abs(volume(B))
norm(B::Scalar) = abs(volume(B))
norm(B::Pseudoscalar) = abs(volume(B))

"""
    sign(B::AbstractBlade)::Int8

Return the sign of a blade relative to its unit basis.
"""
sign(B::Blade)::Int8 = sign(B.volume)
sign(B::Scalar)::Int8 = sign(B.value)
sign(B::Pseudoscalar)::Int8 = sign(B.value)


# --- Core Scalar functions

# Exports
export value

"""
    value(B::Scalar)::Real
    value(B::Pseudoscalar)::Real

Return the value of a scalar.
"""
value(B::Scalar) = B.value
value(B::Pseudoscalar) = B.value


# --- Utility functions

# Exports
export blade_atol

"""
    convert(::Type{S}, B::Blade) where {T<:AbstractFloat, S<:Blade{T}}

Convert Blade to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Blade) where {T<:AbstractFloat, S<:Blade{T}} =
    Blade{T}(B)

"""
    convert(::Type{S}, B::Scalar) where {T<:AbstractFloat, S<:Scalar{T}}

Convert Scalar to have the floating-point precision of type `T`.
"""
convert(::Type{S}, B::Scalar) where {T<:AbstractFloat, S<:Scalar{T}} =
    Scalar{T}(value(B))

#=
 TODO: review numerical error in factorizations to see if a different
       tolerance would be better.
=#
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
