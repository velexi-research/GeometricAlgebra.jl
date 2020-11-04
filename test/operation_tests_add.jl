"""
Unit tests for the +(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

@testset "+(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_plus_C = B + C
    expected_result = Pseudoscalar(test_dim, test_value_1 + test_value_2)
    @test B_plus_C isa Pseudoscalar
    @test B_plus_C == expected_result
end

@testset "+(B, C): B or C isa Scalar" begin
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

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_plus_C = B + C
    expected_result = 1 + test_value_1
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value_1
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    @test B + C == B
    @test C + B == B

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

@testset "+(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    expected_result = 2
    C_plus_B = C + B
    @test B + C isa Scalar
    @test B + C == expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    @test isone(B + C)
    @test isone(C + B)

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_plus_C = B + C
    expected_result = 1 + test_value
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result
end

@testset "+(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    @test iszero(B + C)

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_plus_C = B + C
    expected_result = C
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = C
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result
end
