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

import Base.:(==), Base.:(≈), Base.:(-)


# --- Comparison operators

"""
    ==(B1::AbstractBlade{<:AbstractFloat}, B2::AbstractBlade{<:AbstractFloat})

Return true if B1 and B2 are equal; otherwise, return false.
"""
==(B1::Blade, B2::Blade) =  # TODO: fix basis and orientation comparison
    dim(B1) == dim(B2) && grade(B1) == grade(B2) &&
    value(B1) == value(B2) && basis(B1) == basis(B2)

==(B1::AbstractScalar, B2::AbstractScalar) = (value(B1) == value(B2))

==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (B == x)


"""
    ≈(B1::AbstractBlade{<:AbstractFloat}, B2::AbstractBlade{<:AbstractFloat})

Return true if B1 and B2 are approximatly equal; otherwise, return false.
"""
≈(B1::Blade{T1}, B2::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    sign(B1) == sign(B2) &&
    ≈(norm(B1), norm(B2), atol=atol, rtol=rtol) &&
    true  # TODO: add check that basis represent the same space

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


# --- Unary operations

# Exports
export opposite, reciprocal

"""
    -(B::AbstractBlade), opposite(B::AbstractBlade)

Return the additive inverse of `B`.
"""
opposite(B::Blade{<:AbstractFloat}) = Blade(B, value=-value(B))
opposite(B::Scalar) = Scalar(-value(B))

-(B::AbstractBlade{<:AbstractFloat}) = opposite(B)

"""
    reciprocal(B::AbstractBlade{<:AbstractFloat})

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Blade{<:AbstractFloat}) =
    mod(grade(B), 4) < 2 ?
        Blade(B, value=1 / norm(B)) :
        Blade(B, value=-1 / norm(B))

reciprocal(B::Scalar{<:AbstractFloat}) = Scalar(1 / value(B))
reciprocal(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
reciprocal(B::One{T}) where {T<:AbstractFloat} = One{T}()


# --- Binary operations

# Exports
export ∧, outer

"""
    ∧(B::AbstractBlade, C::AbstractBlade)

    ∧(v::Vector, C::AbstractBlade)
    ∧(B::AbstractBlade, v::Vector)

    ∧(v::Vector, w::Vector)

    ∧(x::Real, C::AbstractBlade)
    ∧(B::AbstractBlade, x::Real)

Return the outer product of a combination of blades and vectors.
"""
∧(B::Blade, C::Blade) =
    Blade(hcat(basis(B), basis(C)), value=value(B) * value(C))

∧(B::AbstractScalar, C::Blade) = Blade(C, value=value(B) * value(C))
∧(B::Blade, C::AbstractScalar) = Blade(B, value=value(B) * value(C))

∧(v::Vector, B::Blade) = Blade(v) ∧ B
∧(B::Blade, v::Vector) = B ∧ Blade(v)
∧(v::Vector, w::Vector) = Blade(v) ∧ Blade(w)

∧(x::Real, B::Blade) = Blade(B, value=x * value(B))
∧(B::Blade, x::Real) = x ∧ B

const outer = ∧

"""
    ⋅(B::Blade, C::Blade)
    ⋅(B::Vector, C::Blade)
    ⋅(B::Blade, C::Vector)
    ⋅(v::Vector, w::Vector)

Return the inner product of a combination of vectors and blades.
"""
# TODO
#⋅(v::Vector{<:AbstractFloat}, w::Vector{<:AbstractFloat}) = 3
#⋅(v::Blade{<:AbstractFloat}, w::Blade{<:AbstractFloat}) = 4
