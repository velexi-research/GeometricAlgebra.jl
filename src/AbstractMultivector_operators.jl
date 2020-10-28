"""
AbstractMultivector_operators.jl defines operators for the AbstractMultivector
type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

# --- Comparison operators

# Note: the comparison operators defined here ensure that arbitrary
#       AbstractMultivector instances can be compared

import Base.:(==), Base.:(≈)

"""
    ==(B::AbstractMultivector, C::AbstractMultivector)

Return true if B and C are equal; otherwise, return false.
"""
==(B::AbstractMultivector, C::AbstractMultivector) = false

"""
    ≈(B::AbstractMultivector, C::AbstractMultivector)

Return true if B and C are approximately equal; otherwise, return false.
"""
≈(B::AbstractMultivector, C::AbstractMultivector) = false
≈(B::AbstractMultivector, C::Real) = false
≈(B::Real, C::AbstractMultivector) = false

# --- Special cases

import Base.:(+), Base.:(-)
import Base.:(*), Base.:(/)

# ------ Operations involving Zero

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
export contractl
contractl(B::Zero, M::AbstractMultivector) = B
contractl(M::AbstractMultivector, B::Zero) = B

# ------ Operations involving One

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
