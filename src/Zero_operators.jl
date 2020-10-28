"""
Zero_oeprators.jl defines operators for the Zero type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Operators from the AbstractMultivector and AbstractBlade interfaces

# ------ Unary operators

-(B::Zero) = B

dual(B::Zero; dim::Union{Integer, Nothing}=nothing) =
    error("The dual of Zero is not well-defined")

# --- Special cases

# Operations involving AbstractMultivectors
+(B::Zero, M::AbstractMultivector) = M
+(M::AbstractMultivector, B::Zero) = M

-(B::Zero, M::AbstractMultivector) = -M
-(M::AbstractMultivector, B::Zero) = M

*(B::Zero, M::AbstractMultivector) = B
*(M::AbstractMultivector, B::Zero) = B

/(B::Zero, M::AbstractMultivector) = B
/(M::AbstractMultivector, B::Zero) = Scalar{typeof(norm(M))}(Inf)

wedge(B::Zero, M::AbstractMultivector) = B
wedge(M::AbstractMultivector, B::Zero) = B

contractl(B::Zero, M::AbstractMultivector) = B
contractl(M::AbstractMultivector, B::Zero) = B

# Operations involving Scalars
+(B::Zero, C::Scalar) = C
+(B::Scalar, C::Zero) = B

-(B::Zero, C::Scalar) = -C
-(B::Scalar, C::Zero) = B

*(B::Zero, C::Scalar) = B
*(B::Scalar, C::Zero) = C

/(B::Zero, C::Scalar) = B
/(B::Scalar, C::Zero) =
    sign(B) > 0 ?
        Scalar{typeof(norm(B))}(Inf) :
        Scalar{typeof(norm(B))}(-Inf)

wedge(B::Zero, C::Scalar) = B
wedge(B::Scalar, C::Zero) = C

contractl(B::Zero, C::Scalar) = B
contractl(B::Scalar, C::Zero) = C

# Operations involving Ones
+(B::Zero, C::One) = C
+(B::One, C::Zero) = B

-(B::Zero, C::One) = -C
-(B::One, C::Zero) = B

*(B::Zero, C::One) = B
*(B::One, C::Zero) = C

/(B::Zero, C::One) = B
/(B::One, C::Zero) = reciprocal(C)

wedge(B::Zero, C::One) = B
wedge(B::One, C::Zero) = C

contractl(B::Zero, C::One) = B
contractl(B::One, C::Zero) = C

# Operations between Zeros
+(B::Zero, C::Zero) = B
-(B::Zero, C::Zero) = B
*(B::Zero, C::Zero) = B
/(B::Zero, C::Zero) = Scalar{typeof(value(B))}(NaN)

wedge(B::Zero, C::Zero) = B
contractl(B::Zero, C::Zero) = B
