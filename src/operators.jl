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
    ==(B1::AbstractBlade, B2::AbstractBlade)

Return true if B1 and B2 are equal; otherwise, return false.
"""
==(B1::Blade, B2::Blade) =
    dim(B1) == dim(B2) && grade(B1) == grade(B2) &&
    volume(B1) == volume(B2) && basis(B1) == basis(B2)

==(B1::Scalar, B2::Scalar) = (value(B1) == value(B2))

==(B::Scalar, x::Real) = (x == value(B))
==(x::Real, B::Scalar) = (B == x)

==(B1::Pseudoscalar, B2::Pseudoscalar) =
    (dim(B1) == dim(B2)) && (value(B1) == value(B2))

"""
    ≈(B1::AbstractBlade, B2::AbstractBlade)

Return true if B1 and B2 are approximatly equal; otherwise, return false.
"""
function ≈(B1::Blade{T1}, B2::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat}
    # Check dim, grade, and norm
    if dim(B1) != dim(B2) || grade(B1) != grade(B2) ||
        ≉(norm(B1), norm(B2), atol=atol, rtol=rtol)

        return false
    end

    # Check that B1 and B2 represent the same space
    projection = LinearAlgebra.det(transpose(basis(B1)) * basis(B2))
    if ≉(abs(projection), 1, atol=atol, rtol=rtol)
        return false
    end

    # Check that B1 and B2 have the same orientation
    return sign(B1) * sign(B2) == sign(projection)
end

≈(B1::Scalar{T1}, B2::Scalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B1), value(B2), atol=atol, rtol=rtol)

≈(B::Scalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)
≈(x::Real, B::Scalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)

≈(B1::Pseudoscalar{T1}, B2::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B1) == dim(B2)) && ≈(value(B1), value(B2), atol=atol, rtol=rtol)

# By default, AbstractBlades are not approximately equal
≈(B1::AbstractBlade, B2::AbstractBlade) = false

# By default, AbstractBlades and Real numbers are not approximately equal
≈(B1::AbstractBlade, B2::Real) = false
≈(B1::Real, B2::AbstractBlade) = false


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

Return the dual `B`.
"""
# TODO: implement
dual(B::Blade) = nothing
dual(B::Scalar) = nothing
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
