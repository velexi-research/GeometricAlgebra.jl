"""
operations.jl defines operations on GeometricAlgebra types

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# Additive operations
include("operations/add.jl")
include("operations/subtract.jl")

# Multiplicative operations
include("operations/multiply.jl")
include("operations/divide.jl")
include("operations/wedge.jl")
include("operations/contractl.jl")

# Geometric operations
include("operations/dual.jl")
include("operations/project.jl")
include("operations/reject.jl")
