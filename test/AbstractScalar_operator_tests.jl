"""
Operator unit tests for concrete subtypes of AbstractScalar.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

using Test
using GeometricAlgebra

# Tests

# --- +(B, C)

@testset "+(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B or C isa Zero

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()
    @test B + C === Zero()

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()
    expected_result = One()
    @test B + C === expected_result
    @test C + B === expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value_2)
    @test B + C === C
    @test C + B === C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value_2
    expected_result = C
    @test B + C isa Scalar
    @test B + C == expected_result
    @test C + B isa Scalar
    @test C + B == expected_result

    # --- B or C isa One

    # B::One, C::One
    B = One()
    C = One()
    @test B + C isa Scalar
    @test B + C == 2

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value_2)
    expected_result = 1 + test_value_2
    @test B + C isa Scalar
    @test B + C == expected_result
    @test C + B isa Scalar
    @test C + B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value_2
    expected_result = 1 + test_value_2
    @test B + C isa Scalar
    @test B + C == expected_result
    @test C + B isa Scalar
    @test C + B == expected_result

    # --- B or C isa Scalar

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)
    @test B + C isa Scalar
    @test B + C == test_value_1 + test_value_2

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2
    expected_result = test_value_1 + test_value_2
    @test B + C isa Scalar
    @test B + C == expected_result
    @test C + B isa Scalar
    @test C + B == expected_result
end
