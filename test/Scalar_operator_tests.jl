"""
Operator unit tests for Scalar type.

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

    # --- Tests

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

    # --- Tests

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

    # --- Tests

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

    # --- Tests

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
