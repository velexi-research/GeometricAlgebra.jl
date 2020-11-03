"""
Unit tests for the dual(x, y) function

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

@testset "dual(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)
        dual_B = dual(B, C)
        expected_result = test_value_1
        @test dual_B isa Scalar
        @test dual_B == expected_result

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) == B
        else
            @test dual(dual_B, C) == -B
        end
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dual(B, C)

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    C = Scalar(test_value_2)

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) === expected_result

        expected_result = mod(dim_B, 4) < 2 ?
            Pseudoscalar(dim_B, test_value_2) :
            Pseudoscalar(dim_B, -test_value_2)
        @test dual(C, B) == expected_result
    end

    # --- B::Pseudoscalar, C::One
    #     B::One, C::Pseudoscalar

    C = One()

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result

        expected_result = mod(dim_B, 4) < 2 ?
            Pseudoscalar(dim_B, 1) :
            Pseudoscalar(dim_B, -1)
        @test dual(C, B) == expected_result
    end

    # --- B::Pseudoscalar, C::Real
    #     B::Real, C::Pseudoscalar

    C = test_value_2

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result
    end
end

@testset "dual(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    B = Scalar(test_value_1)

    # B::Scalar, C::Scalar
    C = Scalar(test_value_2)
    B_dual_C = dual(B, C)
    expected_result = test_value_1
    @test B_dual_C isa Scalar
    @test B_dual_C == expected_result

    # B::Scalar, C::Real
    C = test_value_2
    B_dual_C = dual(B, C)
    expected_result = test_value_1
    @test B_dual_C isa Scalar
    @test B_dual_C == expected_result
end

@testset "dual(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    B = One()

    # B::One, C::One
    C = One()
    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    # B::One, C::Scalar
    # C::Scalar, B::One
    C = Scalar(test_value)

    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    C_dual_B = dual(C, B)
    expected_result = C
    @test C_dual_B == expected_result

    # B::One, C::Real
    # C::Real, C::One
    C = test_value

    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    C_dual_B = dual(C, B)
    expected_result = Scalar(C)
    @test C_dual_B == expected_result
end

@testset "dual(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- dual(B::Zero, C)

    B = Zero()
    expected_error = "The dual of Zero is not well-defined"

    # C::One
    C = One()
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Scalar
    C = Scalar(test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Real
    C = test_value
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Blade
    C = Blade(randn(4, 3))
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Pseudoscalar
    C = Pseudoscalar(5, test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # --- dual(B, C::Zero)

    C = Zero()

    expected_error = "The dual of anything relative to Zero is not well-defined"

    # B::Zero
    B = Zero()
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::One
    B = One()
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Scalar
    B = Scalar(test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Real
    B = test_value
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Blade
    B = Blade(randn(4, 3))
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Pseudoscalar
    B = Pseudoscalar(5, test_value)
    @test_throws ErrorException(expected_error) dual(B, C)
end
