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

# --- Special cases

# Operations between One instances
+(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Scalar{T}(2)
-(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Zero{T}()
*(B::One, C::One) = B
/(B::One, C::One) = B

wedge(B::One, C::One) = B

dot(B::One, C::One) = contractl(B, C)

export contractl
contractl(B::One, C::One) = B

# Operations involving Zero
+(B::One, C::Zero) = B
+(B::Zero, C::One) = C

-(B::One, C::Zero) = B
-(B::Zero, C::One) = -C

*(B::One, C::Zero) = C
*(B::Zero, C::One) = B

/(B::One, C::Zero) = reciprocal(C)
/(B::Zero, C::One) = B

wedge(B::One, C::Zero) = C
wedge(B::Zero, C::One) = B

dot(B::One, C::Zero) = contractl(B, C)
dot(B::Zero, C::One) = contractl(C, B)

export contractl
contractl(B::One, C::Zero) = C
contractl(B::Zero, C::One) = B
