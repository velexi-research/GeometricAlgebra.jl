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

# External packages
using Documenter
using Test
using TestSetExtensions

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Test sets

@testset "Doctests" begin
    doctest(GeometricAlgebra)
end

@testset ExtendedTestSet "Unit tests" begin
    @includetests ARGS
end
