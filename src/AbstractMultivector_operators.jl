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
