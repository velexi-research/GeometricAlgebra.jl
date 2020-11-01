"""
Operator unit tests for One type.

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

@testset "One: -(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        minus_B = -B
        @test minus_B isa Scalar{precision_type}
        @test minus_B == Scalar{precision_type}(-1)
    end
end

@testset "One: reverse(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        reverse_B = reverse(B)
        @test reverse_B === B
        @test B * reverse_B == norm(B)^2
    end
end

@testset "One: dual(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)

            expected_dual = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1) :
                Pseudoscalar{precision_type}(test_dim, -1)
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_dual
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "One: reciprocal(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B === B
    end
end

@testset "One: +(B, C)" begin
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

    expected_result = One()
    @test B + C === expected_result
    @test C + B === expected_result

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

@testset "One: -(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # ------ B::One, C::One

    B = One()
    C = One()

    expected_result = Zero()
    @test B - C === expected_result

    # ------ B::One, C::Zero
    #        B::Zero, C::One

    B = One()
    C = Zero()

    @test B - C === One()

    C_minus_B = C - B
    @test C_minus_B isa Scalar
    @test C_minus_B == -1

    # ------ B::One, C::Real
    #        B::Real, C::One

    # B != C
    B = One()
    C = test_value

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # B == C
    B = One()
    C = 1

    expected_result = Zero()

    B_minus_C = B - C
    @test B_minus_C === expected_result

    C_minus_B = C - B
    @test C_minus_B === expected_result
end

@testset "One: *(B, C)" begin
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

@testset "One: /(B, C)" begin
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

@testset "One: wedge(B, C), B ∧ C" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = One()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = Zero()
    @test C_wedge_B === expected_result
    @test C ∧ B === C_wedge_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_wedge_C = wedge(B, C)
    expected_result = C
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = C
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B
end

@testset "One: contractl(B, C), B < C" begin
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

@testset "One: proj(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_proj_C = proj(B, C)
    expected_result = One()
    @test B_proj_C === expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    expected_result = Zero()

    B_proj_C = proj(B, C)
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    @test C_proj_B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_proj_C = proj(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = C
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result

    # --- B::One, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::One

    # return_blade == true
    B = One()
    C = Vector(rand(10))
    @test proj(B, C) === B
    @test proj(C, B) == zero(B)

    # return_blade == false
    @test proj(C, B, return_blade=false) == 0
    @test proj(B, C, return_blade=false) == value(B)
end

@testset "One: dual(B, C)" begin
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
