"""
AbstractScalar_operators.jl defines the AbstractScalar type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

# --- Operators from the AbstractMultivector and AbstractBlade interfaces

import Base.:(+), Base.:(-)
import Base.:(*), Base.:(/)

+(B::AbstractScalar, C::AbstractScalar) = Scalar(B, value=value(B) + value(C))
-(B::AbstractScalar, C::AbstractScalar) = Scalar(B, value=value(B) - value(C))
*(B::AbstractScalar, C::AbstractScalar) = Scalar(B, value=value(B) * value(C))
/(B::AbstractScalar, C::AbstractScalar) = Scalar(B, value=value(B) / value(C))

# --- Operators between AbstractScalar and Real types

+(B::AbstractScalar, C::Real) = Scalar(B, value=value(B) + C)
+(B::Real, C::AbstractScalar) = C + B

-(B::AbstractScalar, C::Real) = Scalar(B, value=value(B) - C)
-(B::Real, C::AbstractScalar) = Scalar(C, value=B - value(C))

*(B::AbstractScalar, C::Real) = Scalar(B, value=value(B) * C)
*(B::Real, C::AbstractScalar) = C * B

/(B::AbstractScalar, C::Real) = Scalar(B, value=value(B) / C)
/(B::Real, C::AbstractScalar) = Scalar(C, value=B / value(C))

# --- Comparison operators

import Base.:(==), Base.:(≈)

"""
    ==(B::AbstractScalar, C::AbstractScalar)
    ==(B::AbstractScalar, x::Real)
    ==(x::Real, B::AbstractScalar)

Return true if B and C are equal; otherwise, return false.
"""
==(B::AbstractScalar, C::AbstractScalar) = (value(B) == value(C))
==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (value(B) == x)

"""
    ≈(B::AbstractScalar, C::AbstractScalar)
    ≈(B::AbstractScalar, x::Real)
    ≈(x::Real, B::AbstractScalar)

Return true if B and C are approximatly equal; otherwise, return false.
"""
≈(B::AbstractScalar{T1}, C::AbstractScalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B), value(C), atol=atol, rtol=rtol)

≈(B::AbstractScalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)

≈(x::Real, B::AbstractScalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)
