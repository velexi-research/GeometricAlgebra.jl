"""
One_operators.jl defines operators for the One type

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
import LinearAlgebra.dot, LinearAlgebra.:(⋅)

# ------ Unary operators

reciprocal(B::One) = B

# --- Special cases

# Operations involving AbstractMultivector instances
+(M::AbstractMultivector, B::One) = Multivector(vcat([B], blades(M)))
+(B::One, M::AbstractMultivector) = B + M

-(M::AbstractMultivector, B::One) = Multivector(vcat([-B], blades(M)))
-(B::One, M::AbstractMultivector) = -(M - B)

*(M::AbstractMultivector, B::One) = M
*(B::One, M::AbstractMultivector) = M

/(M::AbstractMultivector, B::One) = M
/(B::One, M::AbstractMultivector) = nothing  # TODO

wedge(M::AbstractMultivector, B::One) = M
wedge(B::One, M::AbstractMultivector) = M

contractl(M::AbstractMultivector, B::One) = nothing  # TODO
contractl(B::One, M::AbstractMultivector) = nothing  # TODO

# Operations involving Scalar instances
+(B::Scalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
+(B::One, C::Scalar) = C + B

-(B::Scalar, C::One) = Scalar{typeof(value(B))}(value(B) - 1)
-(B::One, C::Scalar) = -(C - B)

*(B::Scalar, C::One) = B
*(B::One, C::Scalar) = C

/(B::Scalar, C::One) = B
/(B::One, C::Scalar) = reciprocal(C)

wedge(B::Scalar, C::One) = B
wedge(B::One, C::Scalar) = C

contractl(B::Scalar, C::One) = B
contractl(B::One, C::Scalar) = C

# Operations between One instances
+(B::One, C::One) = Scalar{typeof(value(B))}(2)
-(B::One, C::One) = Zero{typeof(value(B))}()
*(B::One, C::One) = B
/(B::One, C::One) = B

wedge(B::One, C::One) = B
contractl(B::One, C::One) = B
