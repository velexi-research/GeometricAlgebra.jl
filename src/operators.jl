"""
The operators.jl submodule defines operations on subtypes of AbstractBlade.

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


# --- Comparison operators

# Imports
import Base.:(==), Base.:(≈)

"""
    ==(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are equal; otherwise, return false.
"""
==(B::Blade, C::Blade) =
    dim(B) == dim(C) && grade(B) == grade(C) &&
    volume(B) == volume(C) && basis(B) == basis(C)

==(B::Scalar, C::Scalar) = (value(B) == value(C))

==(B::Scalar, x::Real) = (x == value(B))
==(x::Real, B::Scalar) = (B == x)

==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

"""
    ≈(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are approximatly equal; otherwise, return false.
"""
function ≈(B::Blade{T1}, C::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat}
    # Check dim, grade, and norm
    if dim(B) != dim(C) || grade(B) != grade(C) ||
        ≉(norm(B), norm(C), atol=atol, rtol=rtol)

        return false
    end

    # Check that B and C represent the same space
    projection = LinearAlgebra.det(transpose(basis(B)) * basis(C))
    if ≉(abs(projection), 1, atol=atol, rtol=rtol)
        return false
    end

    # Check that B and C have the same orientation
    return sign(B) * sign(C) == sign(projection)
end

≈(B::Scalar{T1}, C::Scalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B), value(C), atol=atol, rtol=rtol)

≈(B::Scalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)
≈(x::Real, B::Scalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)

≈(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && ≈(value(B), value(C), atol=atol, rtol=rtol)

# By default, AbstractBlades are not approximately equal
≈(B::AbstractBlade, C::AbstractBlade) = false

# By default, AbstractBlades and Real numbers are not approximately equal
≈(B::AbstractBlade, C::Real) = false
≈(B::Real, C::AbstractBlade) = false


# --- Unary operations

# Imports
import Base.:(-)

# Exports
export dual, reciprocal, reverse

"""
    -(B::AbstractBlade)

Return the additive inverse of `B`.
"""
-(B::Blade{<:AbstractFloat}) = Blade(B, volume=-volume(B))
-(B::Scalar) = Scalar(-value(B))
-(B::Pseudoscalar) = Pseudoscalar(dim(B), -value(B))

"""
    reciprocal(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Blade) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / norm(B)) :
        Blade(B, volume=-1 / norm(B))

reciprocal(B::Scalar) = Scalar(1 / value(B))

reciprocal(B::Pseudoscalar) =
    mod(dim(B), 4) < 2 ?
        Pseudoscalar(dim(B), 1 / value(B)) :
        Pseudoscalar(dim(B), -1 / value(B))

"""
    reverse(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
# TODO: implement


"""
    dual(B::AbstractBlade)

Return the dual `B` (relative to the space that the geometric algebra is
extended from).
"""
dual(B::Blade) = nothing  # TODO: implement
dual(B::Pseudoscalar) = Scalar(value(B))


# --- Binary operations

# Imports
import Base.:(*)

# Exports
export ∧, outer

"""
    ∧(B::Union{AbstractBlade, Vector, Real},
      C::Union{AbstractBlade, Vector, Real})

    outer(B::Union{AbstractBlade, Vector, Real},
          C::Union{AbstractBlade, Vector, Real})

Return the outer product of `B` and `C`.
"""
∧(B::Blade, C::Blade) =
    Blade(hcat(basis(B), basis(C)), volume=volume(B) * volume(C))

∧(B::Scalar, C::Blade) = Blade(C, volume=volume(B) * volume(C))
∧(B::Blade, C::Scalar) = Blade(B, volume=volume(B) * volume(C))

∧(v::Vector, B::Blade) = Blade(v) ∧ B
∧(B::Blade, v::Vector) = B ∧ Blade(v)
∧(v::Vector, w::Vector) = Blade(v) ∧ Blade(w)

∧(x::Real, B::Blade) = Blade(B, volume=x * volume(B))
∧(B::Blade, x::Real) = x ∧ B

outer(x::Union{AbstractBlade, Vector, Real},
      y::Union{AbstractBlade, Vector, Real}) = x ∧ y

"""
    ⋅(B::Union{AbstractBlade, Vector}, C::Union{AbstractBlade, Vector})
    inner(B::Union{AbstractBlade, Vector}, C::Union{AbstractBlade, Vector})

TODO: fill in other function signatures

Return the inner product (left contraction) of `B` and `C`.
"""
# TODO: implement

"""
    *(B::Union{AbstractBlade, Real}, C::Union{AbstractBlade, Real})

TODO: add documentation

Return the geometric product of `B` and `C`.
"""
# TODO: implement
*(x::Real, B::Blade{<:AbstractFloat}) = Blade(B, volume=x * volume(B))
*(B::Blade{<:AbstractFloat}, x::Real) = x * B
*(x::Scalar, B::Blade{<:AbstractFloat}) =
    Blade(B, volume=volume(x) * volume(B))
*(B::Blade{<:AbstractFloat}, x::Scalar) = x * B

"""
    dual(B::AbstractBlade, C::AbstractBlade)

Return the dual `B` relative to `C`.

Notes
-----
* `dual(B, C)` is only defined if (1) `B` and `C` are extended from real
  vector spaces of the same dimension and (2) the subspace represented by `B`
  is contained in subspace represented by `C`.

* The volume of `C` is ignored.
"""
function dual(B::Blade, C::Blade)
    # --- Handle edge cases

    # Check that B and C are extended from the real vector spaces of the same
    # dimension
    if dim(B) != dim(C)
        throw(DimensionMismatch("`dim(B)` not equal to `dim(C)`"))
    end

    # Check that B is contained in C
    projection_coefficients = transpose(basis(C)) * basis(B)
    if LinearAlgebra.norm(projection_coefficients) ≉ grade(B)
        throw(DimensionMismatch("`B` not contained in `C`"))
    end

    # Subspaces represented by B and C are the same
    if grade(B) == grade(C)
        return Scalar(volume(B))
    end

    # --- Compute dual

    permutation = sortperm(sum(abs.(projection_coefficients), dims=2)[:, 1])
    B_ext = hcat(basis(B), basis(C)[:, permutation[1:grade(C) - grade(B)]])
    F = LinearAlgebra.qr(B_ext)
    Blade(Matrix(F.Q)[:, grade(B) + 1:end], volume=volume(B))
end

function dual(B::Pseudoscalar, C::Pseudoscalar)
    # Handle edge cases
    if dim(B) != dim(C)
        throw(DimensionMismatch("`dim(B)` not equal to `dim(C)`"))
    end

    # Compute dual
    Scalar(value(B))
end
