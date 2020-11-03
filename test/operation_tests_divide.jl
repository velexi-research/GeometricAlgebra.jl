"""
Unit tests for the /(x, y) function

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

@testset "/(B, C): B or C isa Pseudoscalar" begin
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

    B_slash_C = B / C
    expected_result = Scalar(test_value_1 / test_value_2)
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    # B::Pseudoscalar, C::Scalar
    # B::Scalar, C::Pseudoscalar
    C = Scalar(test_value_2)

    B = Pseudoscalar(test_dim, test_value_1)
    expected_result = Pseudoscalar(test_dim, test_value_1 / test_value_2)
    @test B / C == expected_result

    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_2 / test_value_1) :
            Pseudoscalar(test_dim, -test_value_2 / test_value_1)

        @test C / B == expected_result
    end

    # B::Pseudoscalar, C::One
    # B::One, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = One()

    expected_result = B
    @test B / C === expected_result

    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, 1/ test_value_1) :
            Pseudoscalar(test_dim, -1/ test_value_1)

        @test C / B == expected_result
    end

    # B::Pseudoscalar, C::Zero
    # B::Zero, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Zero()

    expected_result = sign(B) > 0 ?
        Pseudoscalar(B, value=Inf) :
        Pseudoscalar(B, value=-Inf)
    @test B / C == expected_result

    expected_result = C
    @test C / B === expected_result

    # B::Pseudoscalar, C::Real
    # C::Real, B::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    expected_result = Pseudoscalar(test_dim, test_value_1 / test_value_2)
    @test B / C == expected_result

    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_2 / test_value_1) :
            Pseudoscalar(test_dim, -test_value_2 / test_value_1)

        @test C / B == expected_result
    end
end

@testset "/(B, C): B or C isa Scalar" begin
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

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_slash_C = B / C
    expected_result = test_value_1
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = 1 / test_value_1
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    # B > 0
    B = Scalar(abs(test_value_1))
    B_slash_C = B / C
    expected_result = Inf
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    # B < 0
    B = Scalar(-abs(test_value_1))
    B_slash_C = B / C
    expected_result = -Inf
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = Zero()
    @test C_slash_B === expected_result

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

    # B::Vector, C::Scalar
    B = rand(5)
    C = Scalar(test_value_2)

    expected_result = Blade(B) / C

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == expected_result
end

@testset "/(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_slash_C = B / C
    expected_result = One()
    @test B_slash_C === expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    B_slash_C = B / C
    expected_result = Inf
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = Zero()
    @test C_slash_B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_slash_C = B / C
    expected_result = 1 / test_value
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Vector, C::One
    B = rand(5)
    C = One()

    expected_result = Blade(B)

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == expected_result
end

@testset "/(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(value(B_slash_C))

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_slash_C = B / C
    expected_result = Zero()
    @test B_slash_C === expected_result

    # C > 0
    C = Scalar(abs(test_value))
    C_slash_B = C / B
    expected_result = Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # C < 0
    C = Scalar(-abs(test_value))
    C_slash_B = C / B
    expected_result = -Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Vector, C::Zero
    B = rand(5)
    C = Zero()

    expected_result = Blade(B) / 0

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == expected_result
end
