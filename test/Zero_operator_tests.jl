"""
Operator unit tests for Zero type.

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

@testset "Zero: -(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        @test -B === B
    end
end

@testset "Zero: dual(B)" begin
    B = Zero()
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B)
end

@testset "Zero: reciprocal(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar{precision_type}(Inf)
    end
end

@testset "Zero: +(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    @test B + C === Zero()

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

@testset "Zero: -(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # ------ B::Zero, C::Zero

    B = Zero()
    C = Zero()

    @test B - C === Zero()

    # ------ B::Zero, C::Real
    #        B::Real, C::Zero

    # B != C
    B = Zero()
    C = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value

    C_minus_B = C - B
    @test C_minus_B isa Scalar
    @test C_minus_B == test_value

    # B == C
    B = Zero()
    C = 0

    expected_result = Zero()

    B_minus_C = B - C
    @test B_minus_C === expected_result

    C_minus_B = C - B
    @test C_minus_B === expected_result
end

@testset "Zero: *(B, C)" begin
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
end

@testset "Zero: /(B, C)" begin
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
end

@testset "Zero: wedge(B, C), B ∧ C" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = Zero()
    @test C_wedge_B === expected_result
    @test C ∧ B === C_wedge_B
end

@testset "Zero: contractl(B, C), B < C" begin
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

@testset "Zero: proj(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_proj_C = proj(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_proj_C = proj(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result
end

@testset "Zero: dual(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    B = Zero()

    # C::Zero
    C = Zero()
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)

    # C::One
    C = One()
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)

    # C::Scalar
    C = Scalar(test_value)
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)

    # C::Blade
    C = Blade(randn(4, 3))
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)

    # C::Pseudoscalar
    C = Pseudoscalar(5, test_value)
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)

    # C::Real
    C = test_value
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B, C)
end
