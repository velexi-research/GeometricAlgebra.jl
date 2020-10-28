"""
Scalar_operators.jl defines operators for the Scalar type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

# --- Special cases

# Operations between Scalars
import Base.:(+), Base.:(-), Base.:(*), Base.:(/)

+(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) + value(C))
-(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) - value(C))
*(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) * value(C))
/(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) / value(C))

wedge(B::Scalar, C::Scalar) = B * C
contractl(B::Scalar, C::Scalar) = B * C
proj(B::Scalar, C::Scalar) = B
dual(B::Scalar, C::Scalar) = B

