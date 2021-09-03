"""
Blade.jl defines the Blade type and core methods

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
export Blade

# Functions
export blade_atol

# --- Type definitions

import LinearAlgebra
using LinearAlgebra: det, diag, dot, qr

"""
    struct Blade{T<:AbstractFloat} <: AbstractBlade{T}

Blade represented with the floating-point precision of type `T`. The norm and
orientation of a Blade are encoded by its `volume`. The norm of a Blade is
equal to `abs(volume)` and the orientation of a Blade relative to its `basis`
is equal to `sign(volume)`.

Notes
=====
* The grade of a Blade type is greater than 0 and less than the dimension of
  the space that the blade is embedded in
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

    Note: this inner constructor intended primarily for use by outer
    constructors to enforce constraints.
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
                throw(ArgumentError("`dim` must be nonnegative"))
            end

            if grade < 0
                throw(ArgumentError("`grade` must be nonnegative"))
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
                        message = "columns of `basis` not normalized"
                        throw(ArgumentError(message))
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
            if det(basis) == 0
                return zero(Blade{T})
            elseif det(basis) > 0
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
end

"""
    Blade(vectors::Array{T};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    Blade(vectors::Array{<:Integer};
          volume::Union{Real, Nothing}=nothing,
          atol::Real=blade_atol(Float64))

Construct a Blade from a collection of vectors stored as (1) the columns of a
matrix or (2) a single vector. If the norm of the blade is less than `atol`, a
Scalar representing zero is returned. If the grade of the blade is equal to
the dimension of the space that the blade is embedded in, a Pseudoscalar is
returned.

By default, `vectors` determines the volume of the blade. However, if `volume`
is specified, `vectors` is only used to define the subspace (including
orientation) represented by the blade.

Notes
=====

Orientation
-----------

* _`vectors` contains more than one vector_. When `volume` is positive, the
  orientation of the blade is the same as the orientation of the outer product
  of the columns of `vectors` (taken in order). When `volume` is negative, the
  orientation of the blade is the opposite of the orientation implied by the
  `vectors`.

* _`vectors` contains a single vector `v`_. When `volume` is positive, the
  orientation of the blade is the same as the direction of `v`. When `volume`
  is negative, the orientation of the blade is the opposite of the direction
  of `v`.

Precision
---------

When the precision is not specified, the following rules are applied to set
the precision of the Blade.

* If `vectors` is an Array of floating-point values, the precision of the
  Blade is inferred from the precision of the elements of `vector`.

* If `vectors` is an Array of integers, the precision of the Blade
  defaults to `Float64`.
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
    F = qr(vectors)
    signed_norm = prod(diag(F.R))

    # Return zero if norm is below tolerance
    if abs(signed_norm) < atol
        return zero(Blade{T})
    end

    # Compute orthonormal basis for subspace
    basis = Matrix(F.Q)

    # Compute volume
    volume = isnothing(volume) ? signed_norm : volume * sign(signed_norm)

    # Return new Blade
    Blade{T}(dims[1], dims[2], basis, volume,
             atol=atol, enforce_constraints=false, copy_basis=false)
end

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
    volume = isnothing(volume) ? norm_v : volume

    # Return new Blade
    Blade{T}(length(v), 1, basis, volume, atol=atol,
             enforce_constraints=false, copy_basis=false)
end

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

    Blade{T}(B::Blade;
          volume::Real=volume(B),
          atol::Real=blade_atol(T),
          copy_basis=false) where {T<: AbstractFloat}

Copy constructors. Construct a Blade representing the same space as `B` having
a specified oriented volume relative to `B`. A Scalar representing zero is
returned if the absolute value of `volume` is less than `atol`.

When `copy_basis` is true, the `basis` of the new Blade is a copy of the
`basis` of the original Blade; otherwise, the `basis` of the new Blade is
reference to the `basis` of the original Blade.
"""
Blade{T}(B::Blade;
         volume::Real=volume(B),
         atol::Real=blade_atol(T),
         copy_basis::Bool=false) where {T<:AbstractFloat} =

    Blade{T}(dim(B), grade(B), convert(Matrix{T}, basis(B)), volume,
             atol=atol, enforce_constraints=false, copy_basis=copy_basis)

Blade(B::Blade;
      volume::Real=B.volume,
      atol::Real=blade_atol(typeof(B.volume)),
      copy_basis=false) =
    Blade{typeof(B.volume)}(B, volume=volume, atol=atol, copy_basis=copy_basis)

"""
    Blade(x::Real)
    Blade(x::Scalar)

Convenience constructors that return a Scalar with value `x`.
"""
Blade(x::Real) = Scalar(x)
Blade(x::AbstractScalar) = Scalar(value(x))

# --- Method definitions for AbstractMultivector interface functions

dim(B::Blade) = B.dim

-(B::Blade) = Blade(B, volume=-volume(B), copy_basis=false)

Base.reverse(B::Blade) =
    mod(grade(B), 4) < 2 ?  B : Blade(B, volume=-volume(B), copy_basis=false)

function dual(B::Blade)
    # --- Extend basis(B) to an orthonormal basis for entire space.

    F = qr(basis(B))

    # --- Compute volume of dual

    # Account for orientation of Q relative to orientation of I formed from
    # standard basis
    dual_volume = volume(B) * sign(det(F.Q))

    # Account for orientation of first grade(B) columns of Q relative to
    # orientation of basis(B)
    if prod(diag(F.R)) < 0
        dual_volume = -dual_volume
    end

    # Account for sign of I^{-1} relative to I
    if mod(dim(B), 4) >= 2
        dual_volume = -dual_volume
    end

    # Account for reversals required to eliminate B
    if mod(grade(B), 4) >= 2
        dual_volume = -dual_volume
    end

    Blade{typeof(volume(B))}(dim(B), dim(B) - grade(B),
                             F.Q[:, grade(B) + 1:end], dual_volume,
                             enforce_constraints=false,
                             copy_basis=false)
end

# --- Method definitions for AbstractBlade interface functions

grade(B::Blade) = B.grade

basis(B::Blade) = B.basis

volume(B::Blade) = B.volume

reciprocal(B::Blade) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / volume(B), copy_basis=false) :
        Blade(B, volume=-1 / volume(B), copy_basis=false)

# --- Comparison methods

==(B::AbstractBlade, C::AbstractBlade) =
    dim(B) == dim(C) && grade(B) == grade(C) &&
    volume(B) == volume(C) && basis(B) == basis(C)

function isapprox(B::Blade{T1}, C::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat}
    # Check dim, grade, and norm are equal
    if dim(B) != dim(C) || grade(B) != grade(C) ||
        !isapprox(norm(B), norm(C), atol=atol, rtol=rtol)

        return false
    end

    # Check that B and C represent the same space
    projection = det(transpose(basis(B)) * basis(C))
    if ≉(abs(projection), 1, atol=atol, rtol=rtol)
        return false
    end

    # Check that B and C have the same orientation
    return sign(B) * sign(C) == sign(projection)
end

# --- Utility methods

convert(::Type{T}, B::Blade) where {T<:AbstractFloat} =
    T == typeof(volume(B)) ? B : Blade{T}(B)

#=
 TODO: review numerical error in factorizations to see if a different
       tolerance would be better.
=#
"""
    blade_atol(::Type{T}) where {T<:AbstractFloat}

Return the minimum value of a nonzero blade's norm.
"""
blade_atol(::Type{T}) where {T<:AbstractFloat} = 100*eps(T)
