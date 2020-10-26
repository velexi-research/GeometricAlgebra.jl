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

    B_plus_C = B + C
    expected_result = C
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = C
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

    # --- B or C isa One

    # B::One, C::One
    B = One()
    C = One()

    expected_result = 2
    C_plus_B = C + B
    @test B + C isa Scalar
    @test B + C == expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value_2)

    B_plus_C = B + C
    expected_result = 1 + test_value_2
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value_2
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value_2

    B_plus_C = B + C
    expected_result = 1 + test_value_2
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value_2
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

    # --- B or C isa Scalar

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_plus_C = B + C
    expected_result = test_value_1 + test_value_2
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_plus_C = B + C
    expected_result = test_value_1 + test_value_2
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = test_value_1 + test_value_2
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result
end

@testset "-(B, C)" begin
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
    @test B - C === Zero()

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -1

    @test C - B === One()

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value_2)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value_2

    @test C - B === C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value_2

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value_2

    C_minus_B = C - B
    @test C_minus_B isa Scalar
    @test C_minus_B == test_value_2

    # --- B or C isa One

    # B::One, C::One
    B = One()
    C = One()

    expected_result = Zero()
    @test B - C === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value_2)

    B_minus_C = B - C
    expected_result = 1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value_2 - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value_2

    B_minus_C = B - C
    expected_result = 1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value_2 - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # --- B or C isa Scalar

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value_2 - test_value_1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result
end

@testset "*(B, C)" begin
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

    expected_result = Zero()
    @test B * C === expected_result

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value_2)

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value_2

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # --- B or C isa One

    # B::One, C::One
    B = One()
    C = One()

    expected_result = One()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value_2)

    expected_result = C
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value_2

    B_times_C = B * C
    expected_result = test_value_2
    @test B_times_C isa Scalar
    @test B_times_C == expected_result

    C_times_B = C * B
    expected_result = test_value_2
    @test C_times_B isa Scalar
    @test C_times_B == expected_result

    # --- B or C isa Scalar

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_times_C = B * C
    expected_result = test_value_1 * test_value_2
    @test B_times_C isa Scalar
    @test B_times_C == expected_result

    C_times_B = C * B
    expected_result = test_value_1 * test_value_2
    @test C_times_B isa Scalar
    @test C_times_B == expected_result

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_times_C = B * C
    expected_result = test_value_1 * test_value_2
    @test B_times_C isa Scalar
    @test B_times_C == expected_result

    C_times_B = C * B
    expected_result = test_value_1 * test_value_2
    @test C_times_B isa Scalar
    @test C_times_B == expected_result
end

@testset "/(B, C)" begin
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

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(value(B_slash_C))

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    B_slash_C = B / C
    expected_result = Zero()
    @test B_slash_C === expected_result

    C_slash_B = C / B
    expected_result = Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value_2)

    B_slash_C = B / C
    expected_result = Zero()
    @test B_slash_C === expected_result

    # C > 0
    C = Scalar(abs(test_value_2))
    C_slash_B = C / B
    expected_result = Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # C < 0
    C = Scalar(-abs(test_value_2))
    C_slash_B = C / B
    expected_result = -Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value_2

    B_slash_C = B / C
    expected_result = Zero()
    @test B_slash_C === expected_result

    # C > 0
    C = Scalar(abs(test_value_2))
    C_slash_B = C / B
    expected_result = Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # C < 0
    C = Scalar(-abs(test_value_2))
    C_slash_B = C / B
    expected_result = -Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # --- B or C isa One

    # B::One, C::One
    B = One()
    C = One()

    B_slash_C = B / C
    expected_result = One()
    @test B_slash_C === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value_2)

    B_slash_C = B / C
    expected_result = 1 / test_value_2
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value_2
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value_2

    B_slash_C = B / C
    expected_result = 1 / test_value_2
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value_2
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # --- B or C isa Scalar

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_slash_C = B / C
    expected_result = test_value_1 / test_value_2
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_slash_C = B / C
    expected_result = test_value_1 / test_value_2
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value_2 / test_value_1
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result
end
