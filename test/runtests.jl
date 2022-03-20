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
# TODO: replace local TestSetExtensions.jl with official package after pull request has
#       been accepted

# --- Imports

# Standard library
using Test

# External packages
using Documenter
#using TestSetExtensions
include("../src-external/TestSetExtensions.jl")
ExtendedTestSet = TestSetExtensions.ExtendedTestSet

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Preparations

ENABLE_FAIL_FAST = get(ENV, "JULIA_TEST_FAIL_FAST", "true")
if ENABLE_FAIL_FAST == "true"
    extended_test_set = ExtendedTestSet{Test.FallbackTestSet}
else
    extended_test_set = ExtendedTestSet
end

# --- Test sets

#=
@testset "Doctests" begin
    doctest(GeometricAlgebra)
end
=#

@testset extended_test_set "Unit tests" begin
    @TestSetExtensions.includetests ARGS
end
