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

# --- Operators from the AbstractMultivector and AbstractBlade interfaces

# ------ Binary operators

# dot(B, C)
import LinearAlgebra.dot, LinearAlgebra.:(⋅)
dot(B::Scalar, C::Scalar; left=true) = contractl(B, C)

# contraction operators
import Base.:(<)
export contractl
contractl(B::Scalar, C::Scalar) = B * C
<(B::Scalar, C::Scalar) = contract(B, C)

# --- Special cases

# Operations involving One
dot(B::One, C::Scalar) = C
dot(B::Scalar, C::One) = B
