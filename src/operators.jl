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

import Base.:(==), Base.:(≈)
import Base.:(-)
import Base.:(*)
import LinearAlgebra


# --- Comparison operators

"""
    ==(B1::AbstractBlade, B2::AbstractBlade)

Return true if B1 and B2 are equal; otherwise, return false.
"""
==(B1::Blade, B2::Blade) =
    dim(B1) == dim(B2) && grade(B1) == grade(B2) &&
    volume(B1) == volume(B2) && basis(B1) == basis(B2)

==(B1::AbstractScalar, B2::AbstractScalar) = (value(B1) == value(B2))

==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (B == x)


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
    ≈(x, value(B), rtol=rtol, atol=atol)
≈(x::Real, B::Scalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, rtol=rtol, atol=atol)

≈(B1::Blade, B2::AbstractScalar) = false
≈(B1::AbstractScalar, B2::Blade) = false


# --- Unary operations

# Exports
export reciprocal, reverse

"""
    -(B::AbstractBlade)

Return the additive inverse of `B`.
"""
-(B::Blade{<:AbstractFloat}) = Blade(B, volume=-volume(B))
-(B::Scalar) = Scalar(-value(B))

"""
    reciprocal(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Blade{<:AbstractFloat}) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / norm(B)) :
        Blade(B, volume=-1 / norm(B))

reciprocal(B::Scalar{<:AbstractFloat}) = Scalar(1 / value(B))
reciprocal(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
reciprocal(B::One{T}) where {T<:AbstractFloat} = One{T}()


"""
    reverse(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reverse(B::Blade{<:AbstractFloat}) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / norm(B)) :
        Blade(B, volume=-1 / norm(B))

reverse(B::Scalar{<:AbstractFloat}) = Scalar(1 / value(B))
reverse(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
reverse(B::One{T}) where {T<:AbstractFloat} = One{T}()


# --- Binary operations

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

∧(B::AbstractScalar, C::Blade) = Blade(C, volume=volume(B) * volume(C))
∧(B::Blade, C::AbstractScalar) = Blade(B, volume=volume(B) * volume(C))

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
#⋅(v::Vector{<:AbstractFloat}, w::Vector{<:AbstractFloat}) = 3
#⋅(v::Blade{<:AbstractFloat}, w::Blade{<:AbstractFloat}) = 4
#inner(x::Union{AbstractBlade, Vector, Real},
#      y::Union{AbstractBlade, Vector, Real}) = x ⋅ y

"""
    *(B::Union{AbstractBlade, Real}, C::Union{AbstractBlade, Real})

TODO: add documentation

Return the geometric product of `B` and `C`.
"""
# TODO: implement
*(x::Real, B::Blade{<:AbstractFloat}) = Blade(B, volume=x * volume(B))
*(B::Blade{<:AbstractFloat}, x::Real) = x * B
*(x::AbstractScalar, B::Blade{<:AbstractFloat}) =
    Blade(B, volume=volume(x) * volume(B))
*(B::Blade{<:AbstractFloat}, x::AbstractScalar) = x * B
