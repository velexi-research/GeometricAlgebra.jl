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
# --- Unary operators from the AbstractMultivector and AbstractBlade interfaces

reciprocal(B::One) = B

# --- Binary operators from the AbstractMultivector and AbstractBlade interfaces

# ------ +(B, C)

+(B::One, C::One) = Scalar{typeof(value(B))}(2)

# Operations involving AbstractMultivector
+(M::AbstractMultivector, B::One) = Multivector(vcat([B], blades(M)))
+(B::One, M::AbstractMultivector) = B + M

# Operations involving Scalar
+(B::Scalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
+(B::One, C::Scalar) = C + B

# ------ -(B, C)

-(B::One, C::One) = Zero{typeof(value(B))}()

# Operations involving AbstractMultivector
-(M::AbstractMultivector, B::One) = Multivector(vcat([-B], blades(M)))
-(B::One, M::AbstractMultivector) = -(M - B)

# Operations involving Scalar
-(B::Scalar, C::One) = Scalar{typeof(value(B))}(value(B) - 1)
-(B::One, C::Scalar) = -(C - B)

# ------ *(B, C)

*(B::One, C::One) = B

# Operations involving AbstractMultivector
*(M::AbstractMultivector, B::One) = M
*(B::One, M::AbstractMultivector) = M

# Operations involving Scalar
*(B::Scalar, C::One) = B
*(B::One, C::Scalar) = C

# ------ /(B, C)

/(B::One, C::One) = B

# Operations involving AbstractMultivector
/(M::AbstractMultivector, B::One) = M
/(B::One, M::AbstractMultivector) = nothing  # TODO

# Operations involving Scalar
/(B::Scalar, C::One) = B
/(B::One, C::Scalar) = reciprocal(C)

# ------ wedge(B, C)

wedge(B::One, C::One) = B

# Operations involving AbstractMultivector
wedge(M::AbstractMultivector, B::One) = M
wedge(B::One, M::AbstractMultivector) = M

# Operations involving Scalar
wedge(B::Scalar, C::One) = B
wedge(B::One, C::Scalar) = C

# ------ contractl(B, C)

contractl(B::One, C::One) = B

# Operations involving AbstractMultivector
contractl(M::AbstractMultivector, B::One) = nothing  # TODO
contractl(B::One, M::AbstractMultivector) = nothing  # TODO

# Operations involving Scalar
contractl(B::Scalar, C::One) = B
contractl(B::One, C::Scalar) = C
