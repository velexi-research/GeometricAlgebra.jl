"""
AbstractBlade_operators.jl defines the AbstractBlade type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

# --- Comparison operators

import Base.:(≈)

"""
    ≈(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are approximately equal; otherwise, return false.
"""
≈(B::AbstractBlade, C::AbstractBlade) = false
≈(B::AbstractBlade, C::Real) = false
≈(B::Real, C::AbstractBlade) = false
