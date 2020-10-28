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

# ------ Unary operators

reciprocal(B::One) = B

# ------ Binary operators

+(B::One, M::AbstractMultivector) = Multivector(vcat([B], blades(M)))
+(M::AbstractMultivector, B::One) = B + M

-(M::AbstractMultivector, B::One) = Multivector(vcat([-B], blades(M)))
-(B::One, M::AbstractMultivector) = -(M - B)

*(M::AbstractMultivector, B::One) = M
*(B::One, M::AbstractMultivector) = M

/(M::AbstractMultivector, B::One) = M
/(B::One, M::AbstractMultivector) = nothing  # TODO

# wedge()
export wedge, ∧
wedge(B::One, M::AbstractMultivector) = M
wedge(M::AbstractMultivector, B::One) = M

# dot()
import LinearAlgebra.dot, LinearAlgebra.:(⋅)
dot(B::One, M::AbstractMultivector; left=true) = M
dot(M::AbstractMultivector, B::One; left=true) = M

# --- Operators to remove method ambiguity

# Operators between One instances
+(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Scalar{T}(2)
-(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Zero{T}()
*(B::One, C::One) = B
/(B::One, C::One) = B

wedge(B::One, C::One) = B

dot(B::One, C::One) = contractl(B, C)

import Base.:(<)
export contractl
contractl(B::One, C::One) = B
<(B::One, C::One) = contractl(B, C)

# Operators between One and Zero instances
+(B::One, C::Zero) = B
+(B::Zero, C::One) = C

-(B::One, C::Zero) = B
-(B::Zero, C::One) = -C

*(B::One, C::Zero) = C
*(B::Zero, C::One) = B

/(B::One, C::Zero) = reciprocal(C)
/(B::Zero, C::One) = B

dot(B::One, C::Zero) = contractl(B, C)
dot(B::Zero, C::One) = contractl(C, B)

import Base.:(<)
export contractl
contractl(B::One, C::Zero) = C
contractl(B::Zero, C::One) = B
<(B::One, C::Zero) = contractl(B, C)
<(B::Zero, C::One) = contractl(B, C)

# Operators between One and Scalar instances
dot(B::One, C::Scalar) = C
dot(B::Scalar, C::One) = B

# Operators between One and AbstractScalar instances
+(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
+(B::One, C::AbstractScalar) = C + B

-(B::AbstractScalar, C::One) = B + -C
-(B::One, C::AbstractScalar) = B + -C

*(B::AbstractScalar, C::One) = B
*(B::One, C::AbstractScalar) = C

/(B::AbstractScalar, C::One) = B
/(B::One, C::AbstractScalar) = Scalar{typeof(value(C))}(1 / value(C))
