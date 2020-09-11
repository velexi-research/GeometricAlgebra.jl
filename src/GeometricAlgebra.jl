"""
The GeometricAlgebra.jl module defines geometric algebra types and functions.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
module GeometricAlgebra

# --- Imports

import LinearAlgebra


# --- Submodules

include("blade.jl")


# --- Functions

# Exports
export zero, one


"""
    zero(B::AbstractBlade)

Return Zero (the additive identity).
"""
zero(B::AbstractBlade) = Zero

"""
    one(B::AbstractBlade)

Return One (the multiplicative identity).
"""
one(B::AbstractBlade) = One

end  # End of GeometricAlgebra.jl module
