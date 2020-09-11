"""
Unit tests for the GeometricAlgebra module.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
using Test

# GeometricAlgebra.jl
import GeometricAlgebra.Blade
import GeometricAlgebra.Zero
import GeometricAlgebra.zero
import GeometricAlgebra.One
import GeometricAlgebra.one


# --- Function tests

@testset "zero() tests" begin
    B = Blade([1 2 3])
    @test zero(B) === Zero
end


@testset "one() tests" begin
    B = Blade([1 2 3])
    @test one(B) === One
end
