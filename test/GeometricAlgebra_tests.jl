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
import GeometricAlgebra.One

import GeometricAlgebra.zero
import GeometricAlgebra.one

import GeometricAlgebra.dim
import GeometricAlgebra.grade
import GeometricAlgebra.norm
import GeometricAlgebra.basis
import GeometricAlgebra.inverse

using GeometricAlgebra

# --- Zero type tests

@testset "Zero type tests" begin
    # zero() test
    B = Blade([1 2 3])
    @test zero(B) === Zero

    # Basic blade functions
    B = Zero()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 0
    @test basis(B) === nothing
    @test inverse(B) === NaN
end


# ---  One type tests

@testset "One type tests" begin
    # one() test
    B = Blade([1 2 3])
    @test one(B) === One

    # Basic blade functions
    B = One()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 1
    @test basis(B) === nothing
    @test inverse(B) === One
end
