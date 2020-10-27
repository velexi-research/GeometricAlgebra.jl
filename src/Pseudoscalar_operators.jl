"""
Pseudoscalar_operators.jl defines the Pseudoscalar type

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

# ------ Unary operators

reciprocal(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Pseudoscalar(B, value=1 / value(B)) :
        Pseudoscalar(B, value=-1 / value(B))

# ------ Binary operators

import LinearAlgebra.dot, LinearAlgebra.:(⋅)

function dot(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)

    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))
end

+(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) + value(C))
-(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) - value(C))

*(B::Pseudoscalar, C::Pseudoscalar) = dot(B, C)

/(B::Pseudoscalar, C::Pseudoscalar) =
    Scalar{typeof(value(B))}(value(B) / value(C))

# --- Comparison operators

import Base.:(==), Base.:(≈)

==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

≈(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && ≈(value(B), value(C), atol=atol, rtol=rtol)
