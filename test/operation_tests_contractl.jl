"""
Unit tests for the contractl(B, C) function

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

@testset "contractl(B, C): B or C isa Pseudoscalar" begin
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

    for test_dim = 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_contractl_C isa Scalar
        @test B_contractl_C == expected_result
        @test (B < C) == B_contractl_C
    end

    # --- B::Pseudoscalar, v::Vector
    #     v::Vector, B::Pseudoscalar

    # dim(v) == dim(B)
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        v = Vector(rand(test_dim))

        v_contractl_B = contractl(v, B)
        expected_result = contractl(Blade(v), B)
        @test v_contractl_B ≈ expected_result
        @test v < B == v_contractl_B

        B_contractl_v = contractl(B, v)
        expected_result = zero(B)
        @test B_contractl_v == expected_result
        @test B < v == B_contractl_v
    end

    # dim(v) != dim(B)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    v = test_vector

    @test_throws DimensionMismatch contractl(v, B)
    @test_throws DimensionMismatch v < B

    @test_throws DimensionMismatch contractl(B, v)
    @test_throws DimensionMismatch B < v
end

@testset "contractl(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    # --- B::Scalar, C::One
    #     B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_contractl_C = contractl(B, C)
    expected_result = B
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = B
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B

    # --- B::Scalar, C::Zero
    #     B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B

    # --- B::Scalar, C::Real
    #     B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B

    # --- B::Scalar, v::Vector
    #     v::Vector, B::Scalar

    B = Scalar(test_value_1)
    v = test_vector

    v_contractl_B = contractl(v, B)
    expected_result = zero(B)
    @test v_contractl_B == expected_result
    @test (v < B) == v_contractl_B

    B_contractl_v = contractl(B, v)
    expected_result = Blade(v, volume=norm(v) * test_value)
    @test B_contractl_v == expected_result
    @test (B < v) == B_contractl_v
end

@testset "contractl(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_contractl_C = contractl(B, C)
    expected_result = One()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_contractl_C = contractl(B, C)
    expected_result = C
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = C
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B
end

@testset "contractl(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B
end
