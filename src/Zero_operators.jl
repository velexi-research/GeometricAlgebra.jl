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

# --- Special cases

# Operations between Zero instances
+(B::Zero, C::Zero) = B
-(B::Zero, C::Zero) = B
*(B::Zero, C::Zero) = B
/(B::Zero{T}, C::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(NaN)

wedge(B::Zero, C::Zero) = B

dot(B::Zero, C::Zero) = contractl(B, C)
contractl(B::Zero, C::Zero) = B
<(B::Zero, C::Zero) = contractl(B, C)
