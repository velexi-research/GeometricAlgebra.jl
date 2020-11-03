"""
Unit tests for the *(x, y) function

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

@testset "*(B, C): B or C isa Blade" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test B * C ≈ expected_result
    @test C * B ≈ expected_result

    # --- B::Real, C::Blade
    #     B::Blade, C::Real

    # Preparations
    x = test_value_1
    B = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x * B ≈ expected_result
    @test B * x ≈ expected_result
end

@testset "*(B, C): B or C isa Pseudoscalar" begin
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
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_times_C = B * C
        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_times_C isa Scalar
        @test B_times_C == expected_result
    end

    # B::Pseudoscalar, C::Scalar
    # B::Scalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B * C == expected_result
    @test C * B == expected_result

    # B::Pseudoscalar, C::One
    # B::One, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = One()

    expected_result = B
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Pseudoscalar, C::Zero
    # B::Zero, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Zero()

    expected_result = C
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Pseudoscalar, C::Real
    # C::Real, B::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B * C == expected_result
    @test C * B == expected_result
end

@testset "*(B, C): B or C isa Scalar" begin
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

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    expected_result = B
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

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

    # B::Vector, C::Scalar
    # B::Scalar, B::Vector
    B = rand(5)
    C = Scalar(test_value_2)

    expected_result = test_value_2 * Blade(B)

    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C ≈ expected_result

    C_times_B = C * B
    @test C_times_B isa Blade
    @test C_times_B ≈ expected_result
end

@testset "*(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    expected_result = One()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_times_C = B * C
    expected_result = test_value
    @test B_times_C isa Scalar
    @test B_times_C == expected_result

    C_times_B = C * B
    expected_result = test_value
    @test C_times_B isa Scalar
    @test C_times_B == expected_result

    # B::Vector, C::One
    # B::One, B::Vector
    B = rand(5)
    C = One()

    expected_result = Blade(B)

    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C == expected_result

    C_times_B = C * B
    @test C_times_B isa Blade
    @test C_times_B == expected_result
end

@testset "*(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    expected_result = Zero()
    @test B * C === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Vector, C::Zero
    # B::Zero, B::Vector
    B = rand(5)
    C = Zero()

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result
end
