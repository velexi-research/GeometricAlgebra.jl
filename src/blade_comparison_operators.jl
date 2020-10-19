"""
The blade_comparison_operators.jl submodule defines operations on subtypes of
AbstractBlade.

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
# Scalar comparison
==(B::Scalar, C::Scalar) = (value(B) == value(C))
==(B::Scalar, x::Real) = (x == value(B))
==(x::Real, B::Scalar) = (B == x)

# Blade comparison
==(B::Blade, C::Blade) =
    dim(B) == dim(C) && grade(B) == grade(C) &&
    volume(B) == volume(C) && basis(B) == basis(C)

# Pseudoscalar comparison
==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

"""
    ≈(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are approximatly equal; otherwise, return false.
"""
# Scalar comparison
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

# Blade comparison
function ≈(B::Blade{T1}, C::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat}
    # Check dim, grade, and norm are equal
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

# Pseudoscalar comparison
≈(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && ≈(value(B), value(C), atol=atol, rtol=rtol)

# Default AbstractBlades comparisons
≈(B::AbstractBlade, C::AbstractBlade) = false
≈(B::AbstractBlade, C::Real) = false
≈(B::Real, C::AbstractBlade) = false
