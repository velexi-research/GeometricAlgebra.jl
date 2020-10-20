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
export Scalar, Blade, Pseudoscalar

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
    Construct a Scalar having the specified `value`.
    """
    Scalar{T}(value::Real) where {T<:AbstractFloat} = new(T(value))
end

"""
    Scalar(value::Real)

Construct a Scalar having the specified `value`.

When the precision is not specified, the following rules are applied to set
the precision of the Scalar.

* If `value` is a floating-point value, the precision of the Scalar is inferred
  from the precision of `value`.

* If `value` is an integer, the precision of the Scalar defaults to `Float64`.
"""
Scalar(value::AbstractFloat) = Scalar{typeof(value)}(value)
Scalar(value::Integer) = Scalar(Float64(value))

"""
    Scalar(B::Scalar; value::Real=value(B))

Copy constructor. Construct a Scalar with the same precision as `B` having the
specified `value`.
"""
Scalar(B::Scalar; value::Real=value(B)) = Scalar{typeof(B.value)}(value)


# Blade
"""
    struct Blade{T<:AbstractFloat} <: AbstractBlade{T}

Blade represented with the floating-point precision of type `T`. The norm and
orientation of a Blade are encoded by its `volume`. The norm of a Blade is
equal to `abs(volume)` and the orientation of a Blade relative to its `basis`
is equal to `sign(volume)`.

Notes
-----
* The grade of a Blade type is greater than 0 and less than the dimension of
  the space  that the blade is embedded in
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
    Construct a Blade for a geometric algebra in `dim` dimensions having the
    specified data field values. If the absolute value of `volume` is less than
    `atol`, a Scalar representing zero is returned. If the `dim` and `grade`
    are equal, a Pseudoscalar is returned.

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
            if dim < 0
                error("`dim` must be nonnegative")
            end

            if grade < 0
                error("`grade` must be nonnegative")
            end

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

            if dim > 0 && grade > 0
                for i in 1:grade
                    if LinearAlgebra.norm(basis[:, i]) ≉ 1
                        error("columns of `basis` not normalized")
                    end
                end
            end
        end

        # --- Handle edge cases

        # `volume` is effectively zero
        if abs(volume) < atol
            return zero(Blade{T})
        end

        # Return zero if dim or grade is equal to 0
        if dim == 0 || grade == 0
            return zero(Blade{T})
        end

        # Return a Pseudoscalar if the grade of the blade is equal to the
        # dimension of the embedding space.
        if grade == dim
            if LinearAlgebra.det(basis) == 0
                return zero(Blade{T})
            elseif LinearAlgebra.det(basis) > 0
                return Pseudoscalar{T}(dim, volume)
            else
                return Pseudoscalar{T}(dim, -volume)
            end
        end

        # Return new Blade
        copy_basis ?
            new(dim, grade, copy(basis), volume) :
            new(dim, grade, basis, volume)
    end

    """
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
    function Blade{T}(vectors::Matrix{<:Real};
                      volume::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `volume` is effectively zero
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
        end

        # number of vectors == 0
        dims = size(vectors)
        if dims[1] == 0 || dims[2] == 0
            return zero(Blade{T})
        end

        # number of vectors > dimension of column space
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

        # Convert `vectors` to have type `T`
        vectors = convert(Matrix{T}, vectors)

        # Compute QR factorization and signed norm
        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))

        # Return zero if norm is below tolerance
        if abs(signed_norm) < atol
            return zero(Blade{T})
        end

        # Compute orthonormal basis for subspace
        basis = Matrix(F.Q)

        # Compute volume
        volume = (volume == nothing) ? signed_norm : volume * sign(signed_norm)

        # Return new Blade
        Blade{T}(dims[1], dims[2], basis, volume,
                 atol=atol, enforce_constraints=false, copy_basis=false)
    end

    """
    Construct a Blade from a single vector. A Scalar representing zero is
    returned when the norm of the blade is less than `atol`.

    By default, `v` determines the volume of the blade. However, if `volume`
    is specified, `v` is only used to define the subspace represented by the
    blade.

    When `volume` is positive, the orientation of the blade is the same as the
    direction of `v`. When `volume` is negative, the orientation of the blade
    is the opposite of the direction of `v`.
    """
    function Blade{T}(v::Vector{<:Real};
                      volume::Union{Real, Nothing}=nothing,
                      atol::Real=blade_atol(T)) where {T<:AbstractFloat}

        # --- Handle edge cases

        # `volume` is effectively zero
        if volume != nothing && abs(volume) < atol
            return zero(Blade{T})
        end

        # length(v) == 0
        if length(v) == 0
            return zero(Blade{T})
        end

        # --- Construct Blade

        # Convert `v` to have type `T` and compute norm
        v = convert(Vector{T}, v)
        norm_v = LinearAlgebra.norm(v)

        # Return zero if norm is below tolerance
        if abs(norm_v) < atol
            return zero(Blade{T})
        end

        # Compute basis as a column vector
        basis = reshape(v, length(v), 1) / norm_v

        # Compute volume
        volume = (volume == nothing) ? norm_v : volume

        # Return new Blade
        Blade{T}(length(v), 1, basis, volume, atol=atol,
                 enforce_constraints=false, copy_basis=false)
    end

    """
    Construct a Blade representing the same space as `B` having a specified
    oriented `volume` relative to `B`. A Scalar representing zero is returned
    if the absolute value of `volume` is less than `atol`,

    When `copy_basis` is true, the `basis` of the new Blade is a copy
    of the `basis` of the original Blade; otherwise, the `basis` of the new
    Blade is reference to the `basis` of the original Blade.
    """
    Blade{T}(B::Blade;
             volume::Real=volume(B),
             atol::Real=blade_atol(T),
             copy_basis::Bool=false) where {T<:AbstractFloat} =

        Blade{T}(dim(B), grade(B), convert(Matrix{T}, basis(B)), volume,
                 atol=atol, enforce_constraints=false, copy_basis=copy_basis)
end

"""
    Blade(vectors::Array{T};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(Float64))

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

Blade(vectors::Array{<:Integer};
      volume::Union{Real, Nothing}=nothing,
      atol::Real=blade_atol(Float64)) =
    Blade(convert(Array{Float64}, vectors), volume=volume, atol=atol)

"""
    Blade(B::Blade;
          volume::Real=volume(B),
          atol::Real=blade_atol(typeof(volume(B))),
          copy_basis=false)

Copy constructor. Construct a Blade representing the same space as `B` having
a specified oriented volume relative to `B`. A Scalar representing zero is
returned if the absolute value of `volume` is less than `atol`.

When `copy_basis` is true, the `basis` of the new Blade is a copy of the
`basis` of the original Blade; otherwise, the `basis` of the new Blade is
reference to the `basis` of the original Blade.
"""
Blade(B::Blade;
      volume::Real=B.volume,
      atol::Real=blade_atol(typeof(B.volume)),
      copy_basis=false) =
    Blade{typeof(B.volume)}(B, volume=volume, atol=atol, copy_basis=copy_basis)

"""
    Blade(x::Real)

Convenience constructor that returns a Scalar with value `x`.
"""
Blade(x::Real) = Scalar(x)


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
    Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having
    the specified `value`.
    """
    Pseudoscalar{T}(dim::Integer, value::Real) where {T<:AbstractFloat} =
        dim <= 0 ? error("`dim` must be positive") :
            value == 0 ? zero(Pseudoscalar{T}) : new(dim, T(value))
end

"""
    Pseudoscalar(dim::Integer, value::Real)

Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having the
specified `value`.

When the precision is not specified, the following rules are applied to set
the precision of the Pseudoscalar.

* If `value` is a floating-point value, the precision of the Pseudoscalar is
  inferred from the precision of `value`.

* If `value` is an integer, the precision of the Pseudoscalar defaults to
  `Float64`.
"""
Pseudoscalar(dim::Integer, value::AbstractFloat) =
    Pseudoscalar{typeof(value)}(dim, value)

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))


"""
    Pseudoscalar(B::Pseudoscalar{T};
                 value::Real=value(B)) where {T<:AbstractFloat}

Copy constructor. Construct a Pseudoscalar representing the same space as
`B` having the specified `value`.
"""
Pseudoscalar(B::Pseudoscalar; value::Real=value(B)) =
    Pseudoscalar{typeof(B.value)}(dim(B), value)


# --- Special number functions

# Imports
import Base.zero, Base.one

"""
    zero(M::AbstractMultivector)
    zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(M::AbstractMultivector) = Scalar{typeof(norm(M))}(0)
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Scalar{T}(0)

"""
    one(M::AbstractMultivector)
    one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(M::AbstractMultivector) = Scalar{typeof(norm(M))}(1)
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Scalar{T}(1)


# --- AbstractBlade interface functions: Scalar/Blade/Pseudoscalar

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


# --- Scalar/Pseudoscalar functions

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
    T == typeof(value(B)) ? B : Scalar{T}(value(B))

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
    T == typeof(value(B)) ? B : Pseudoscalar{T}(dim(B), value(B))

#=
 TODO: review numerical error in factorizations to see if a different
       tolerance would be better.
=#
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
