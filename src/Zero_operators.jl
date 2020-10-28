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

import Base.:(+), Base.:(-)
import Base.:(*), Base.:(/)

# ------ Unary operators

-(B::Zero) = B
dual(B::Zero; dim::Union{Integer, Nothing}=nothing) =
    error("The dual of Zero is not well-defined")

# ------ Binary operators

+(B::Zero, M::AbstractMultivector) = M
+(M::AbstractMultivector, B::Zero) = M

-(B::Zero, M::AbstractMultivector) = -M
-(M::AbstractMultivector, B::Zero) = M

*(B::Zero, M::AbstractMultivector) = B
*(M::AbstractMultivector, B::Zero) = B

/(B::Zero, M::AbstractMultivector) = B
/(M::AbstractMultivector, B::Zero) = Scalar{typeof(norm(M))}(Inf)

# wedge()
export wedge, ∧
wedge(B::Zero, M::AbstractMultivector) = B
wedge(M::AbstractMultivector, B::Zero) = B

# dot()
import LinearAlgebra.dot, LinearAlgebra.:(⋅)
dot(B::Zero, M::AbstractMultivector; left=true) = B
dot(M::AbstractMultivector, B::Zero; left=true) = B

# contraction operators
import Base.:(<)
export contractl
contractl(B::Zero, M::AbstractMultivector) = B
contractl(M::AbstractMultivector, B::Zero) = B
<(M::AbstractMultivector, B::Zero) = contractl(B, M)
<(B::Zero, M::AbstractMultivector) = contractl(M, B)

# --- Operators to remove method ambiguity

# Operators between Zero instances
+(B::Zero, C::Zero) = B
-(B::Zero, C::Zero) = B
*(B::Zero, C::Zero) = B
/(B::Zero{T}, C::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(NaN)

wedge(B::Zero, C::Zero) = B

dot(B::Zero, C::Zero) = contractl(B, C)
contractl(B::Zero, C::Zero) = B
<(B::Zero, C::Zero) = contractl(B, C)

# Operators between Zero and AbstractScalar instances
+(B::Zero, C::AbstractScalar) = C
+(B::AbstractScalar, C::Zero) = B

-(B::Zero, C::AbstractScalar) = -C
-(B::AbstractScalar, C::Zero) = B

*(B::Zero, C::AbstractScalar) = B
*(B::AbstractScalar, C::Zero) = C

/(B::Zero, C::AbstractScalar) = B
/(B::AbstractScalar, C::Zero) =
    sign(B) > 0 ?
        Scalar{typeof(norm(B))}(Inf) :
        Scalar{typeof(norm(B))}(-Inf)

dot(B::Zero, C::AbstractScalar) = B
dot(B::AbstractScalar, C::Zero) = C
