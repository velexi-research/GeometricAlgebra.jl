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

include("AbstractMultivector.jl")

include("AbstractBlade.jl")
include("AbstractBlade_operators.jl")

include("AbstractScalar.jl")
include("AbstractScalar_operators.jl")

include("Zero.jl")
include("Zero_operators.jl")

include("One.jl")
include("One_operators.jl")

include("Scalar.jl")

include("Pseudoscalar.jl")

include("Blade.jl")
include("Blade_comparison_operators.jl")
include("Blade_operators.jl")
include("Blade_operators_mgs.jl")

include("Multivector.jl")
include("Multivector_operators.jl")


end  # End of GeometricAlgebra.jl module
