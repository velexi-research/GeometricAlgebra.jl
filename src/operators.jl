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
    B1.dim == B2.dim && B1.grade == B2.grade &&
    B1.norm == B2.norm &&
    B1.basis == B2.basis && B1.sign == B2.sign

==(B1::Scalar, B2::Scalar) = (B1.value == B2.value)

==(B::Scalar, x::Real) = (x == B.value)
==(x::Real, B::Scalar) = (B == x)

==(B1::Scalar, B2::One) = (B1.value == 1)
==(B1::One, B2::Scalar) = (B2 == B1)

==(B1::Zero, B2::Zero) = true
==(B::Zero, x::Real) = (x == 0)
==(x::Real, B::Zero) = (B == x)

==(B1::One, B2::One) = true
==(B::One, x::Real) = (x == 1)
==(x::Real, B::One) = (B == x)

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


# --- Unary operations

# Exports
export negation, reciprocal

"""
    -(B::AbstractBlade), negation(B::AbstractBlade)

Return the additive inverse of `B`.
"""
negation(B::Blade{<:AbstractFloat}) = Blade(B, norm=norm(B), sign=-B.sign)
negation(B::Scalar) = Scalar(-B.value)

-(B::AbstractBlade{<:AbstractFloat}) = negation(B)

"""
    reciprocal(B::AbstractBlade{<:AbstractFloat})

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Blade{<:AbstractFloat}) =
    mod(grade(B), 4) < 2 ?
        Blade(B, norm=1 / norm(B), sign=B.sign) :
        Blade(B, norm=1 / norm(B), sign=-B.sign)

reciprocal(B::Scalar{<:AbstractFloat}) = Scalar(1 / B.value)
reciprocal(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)
reciprocal(B::One{T}) where {T<:AbstractFloat} = One{T}()


# --- Binary operations

# Exports
export ∧, outer

"""
    ∧(B::Blade, C::Blade)
    ∧(B::Blade, C::Vector)
    ∧(B::Vector, C::Blade)
    ∧(v::Vector, w::Vector)

Return the outer product of a combination of vectors and blades.
"""
function ∧(B::Blade{T}, C::Blade{T}) where {T<:AbstractFloat}
    outer_product = Blade{T}(hcat(basis(B), basis(C)))
    if outer_product == 0
        return Zero(T)
    else
        return Blade{T}(outer_product,
                        norm=B.norm * C.norm,
                        sign=outer_product.sign * B.sign * C.sign)
    end
end

∧(B::Blade{T}, C::Vector{T}) where {T<:AbstractFloat} = B ∧ Blade(C)
∧(B::Vector{T}, C::Blade{T}) where {T<:AbstractFloat} = Blade(B) ∧ C
∧(v::Vector{<:AbstractFloat}, w::Vector{<:AbstractFloat}) = Blade(v) ∧ Blade(w)
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