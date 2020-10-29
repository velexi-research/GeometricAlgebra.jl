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
# --- Unary operators from the AbstractMultivector and AbstractBlade interfaces

-(B::Zero) = B

dual(B::Zero; dim::Union{Integer, Nothing}=nothing) =
    error("The dual of Zero is not well-defined")

# --- Binary operators from the AbstractMultivector and AbstractBlade interfaces

# ------ +(B, C)

+(B::Zero, C::Zero) = B

# Operations involving AbstractMultivectors
+(B::Zero, M::AbstractMultivector) = M
+(M::AbstractMultivector, B::Zero) = M

# Operations involving Scalars
+(B::Zero, C::Scalar) = C
+(B::Scalar, C::Zero) = B

# Operations involving Ones
+(B::Zero, C::One) = C
+(B::One, C::Zero) = B

# ------ -(B, C)

-(B::Zero, C::Zero) = B

# Operations involving AbstractMultivectors
-(B::Zero, M::AbstractMultivector) = -M
-(M::AbstractMultivector, B::Zero) = M

# Operations involving Scalars
-(B::Zero, C::Scalar) = -C
-(B::Scalar, C::Zero) = B

# Operations involving Ones
-(B::Zero, C::One) = -C
-(B::One, C::Zero) = B

# ------ *(B, C)

*(B::Zero, C::Zero) = B

# Operations involving AbstractMultivectors
*(B::Zero, M::AbstractMultivector) = B
*(M::AbstractMultivector, B::Zero) = B

# Operations involving Scalars
*(B::Zero, C::Scalar) = B
*(B::Scalar, C::Zero) = C

# Operations involving Ones
*(B::Zero, C::One) = B
*(B::One, C::Zero) = C

# ------ /(B, C)

/(B::Zero, C::Zero) = Scalar{typeof(value(B))}(NaN)

# Operations involving AbstractMultivectors
/(B::Zero, M::AbstractMultivector) = B
/(M::AbstractMultivector, B::Zero) = Scalar{typeof(norm(M))}(Inf)

# Operations involving Scalars
/(B::Zero, C::Scalar) = B
/(B::Scalar, C::Zero) =
    sign(B) > 0 ?
        Scalar{typeof(norm(B))}(Inf) :
        Scalar{typeof(norm(B))}(-Inf)

# Operations involving Ones
/(B::Zero, C::One) = B
/(B::One, C::Zero) = reciprocal(C)

# ------ wedge(B, C)

wedge(B::Zero, C::Zero) = B

# Operations involving AbstractMultivectors
wedge(B::Zero, M::AbstractMultivector) = B
wedge(M::AbstractMultivector, B::Zero) = B

# Operations involving Scalars
wedge(B::Zero, C::Scalar) = B
wedge(B::Scalar, C::Zero) = C

# Operations involving Ones
wedge(B::Zero, C::One) = B
wedge(B::One, C::Zero) = C

# ------ contractl(B, C)

contractl(B::Zero, C::Zero) = B

# Operations involving AbstractMultivectors
contractl(B::Zero, M::AbstractMultivector) = B
contractl(M::AbstractMultivector, B::Zero) = B

# Operations involving Scalars
contractl(B::Zero, C::Scalar) = B
contractl(B::Scalar, C::Zero) = C

# Operations involving Ones
contractl(B::Zero, C::One) = B
contractl(B::One, C::Zero) = C
