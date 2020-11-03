"""
Unit tests for GeometricAlgebra.jl package.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

using Test, TestSetExtensions


# --- Test sets

@testset ExtendedTestSet "All the tests" begin
    @includetests ARGS
end
