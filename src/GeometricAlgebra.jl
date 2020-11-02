"""
The GeometricAlgebra.jl module defines geometric algebra types and functions.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
module GeometricAlgebra

# --- Submodules

# Types
include("AbstractMultivector.jl")
include("AbstractBlade.jl")
include("AbstractScalar.jl")
include("Zero.jl")
include("One.jl")
include("Scalar.jl")
include("Blade.jl")
include("Pseudoscalar.jl")
include("Multivector.jl")

# Methods
include("operations.jl")

end  # End of GeometricAlgebra.jl module
