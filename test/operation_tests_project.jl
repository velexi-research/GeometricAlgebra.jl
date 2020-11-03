"""
Unit tests for the project(x, y) function

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

@testset "project(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # --- B::Pseudoscalar, C::Pseudoscalar

    for test_dim in 5:8
        # dim(B) == dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_proj_C = project(B, C)
        expected_result = B
        @test B_proj_C == expected_result

        # dim(B) != dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim + 1, test_value_2)
        @test_throws DimensionMismatch project(B, C)
    end

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = test_value_1
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == expected_result

    expected_result = zero(C)
    C_proj_B = project(C, B)
    @test C_proj_B === expected_result
end

@testset "project(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 15

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vector = rand(test_dim)

    # --- B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_proj_C = project(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    # --- B::Scalar, C::One
    #     B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_proj_C = project(B, C)
    expected_result = B
    @test B_proj_C == expected_result

    C_proj_B = project(C, B)
    expected_result = C
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Zero
    #     B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Real
    #     B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_proj_C = project(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    C_proj_B = project(C, B)
    expected_result = test_value_2
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result

    # --- B::Scalar, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::Scalar

    # return_blade == true
    B = Scalar(test_value_1)
    C = test_vector
    @test project(B, C) == B
    @test project(C, B) === zero(B)

    # return_blade == false
    @test project(C, B, return_blade=false) == 0
    @test project(B, C, return_blade=false) == value(B)
end

@testset "project(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_proj_C = project(B, C)
    expected_result = One()
    @test B_proj_C === expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    expected_result = Zero()

    B_proj_C = project(B, C)
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    @test C_proj_B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_proj_C = project(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = C
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result

    # --- B::One, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::One

    # return_blade == true
    B = One()
    C = Vector(rand(10))
    @test project(B, C) === B
    @test project(C, B) == zero(B)

    # return_blade == false
    @test project(C, B, return_blade=false) == 0
    @test project(B, C, return_blade=false) == value(B)
end

@testset "project(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # B::Zero, C::Vector
    # B::Vector, C::Zero
    B = Zero()
    C = Vector(rand(10))

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result
end
