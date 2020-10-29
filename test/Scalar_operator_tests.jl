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
# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

@testset "Scalar: -(B)" begin  # from AbstractMultivector interface
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        minus_B = -B
        @test minus_B isa Scalar{precision_type}
        @test minus_B == Scalar{precision_type}(-test_value)
    end
end

@testset "Scalar: dual(B)" begin  # from AbstractMultivector interface
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)

            expected_dual = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(test_value)) :
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(-test_value))
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_dual
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "Scalar: reciprocal(B)" begin  # from AbstractBlade interface
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B ==
            Scalar{precision_type}(1 / precision_type(test_value))
    end
end

@testset "Scalar: +(B, C)" begin
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

    @test B + C === B
    @test C + B === B

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

@testset "Scalar: -(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # ------ B::Scalar, C::Scalar

    # B != C
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = Scalar(test_value_1)
    C = Scalar(test_value_1)

    B_minus_C = B - C
    expected_result = Zero()
    @test B_minus_C === expected_result

    # ------ B::Scalar, C::One
    #        C::One, B::Scalar

    B = Scalar(test_value_1)
    C = One()

    B_minus_C = B - C
    expected_result = test_value_1 - 1
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = 1 - test_value_1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # ------ B::Scalar, C::Zero
    #        B::Zero, C::Scalar

    B = Scalar(test_value_1)
    C = Zero()

    @test B - C === B

    C_minus_B = C - B
    @test C_minus_B isa Scalar
    @test C_minus_B == -test_value_1

    # ------ B::Scalar, C::Real
    #        B::Real, C::Scalar

    # B != C
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

    # B == C
    B = Scalar(test_value_1)
    C = test_value_1

    B_minus_C = B - C
    expected_result = Zero()
    @test B_minus_C === expected_result

    C_minus_B = C - B
    expected_result = Zero()
    @test C_minus_B === expected_result
end

@testset "Scalar: *(B, C)" begin
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
end

@testset "Scalar: /(B, C)" begin
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
end

@testset "Scalar: wedge(B, C), B ∧ C" begin
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

    B_wedge_C = wedge(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C === B_wedge_C

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = B
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = B
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B === C_wedge_B

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = Zero()
    @test C_wedge_B === expected_result
    @test C ∧ B === C_wedge_B

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_wedge_C = wedge(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B === C_wedge_B
end

@testset "Scalar: contractl(B, C), B < C" begin
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

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) === B_contractl_C

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_contractl_C = contractl(B, C)
    expected_result = B
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = B
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) === C_contractl_B

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
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

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) === C_contractl_B
end

@testset "Scalar: proj(B, C)" begin
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

    B_proj_C = proj(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_proj_C = proj(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = C
    @test C_proj_B === expected_result

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_proj_C = proj(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_proj_C = proj(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    C_proj_B = proj(C, B)
    expected_result = test_value_2
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result
end

@testset "Scalar: dual(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    B = Scalar(test_value_1)

    # C::Scalar
    C = Scalar(test_value_2)
    B_dual_C = dual(B, C)
    expected_result = test_value_1
    @test B_dual_C isa Scalar
    @test B_dual_C == expected_result

    # C::Blade
    test_dim = 10
    for test_grade in 5:8
        C = Blade(randn(test_dim, test_grade))
        B_dual_C = dual(B, C)
        expected_result = mod(test_grade, 4) < 2 ?
            Blade(C, volume=test_value_1) :
            Blade(C, volume=-test_value_1)
        @test B_dual_C isa Blade
        @test B_dual_C == expected_result
    end

    # C::Pseudoscalar
    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_2)
        B_dual_C = dual(B, C)
        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_1) :
            Pseudoscalar(test_dim, -test_value_1)
        @test B_dual_C isa Pseudoscalar
        @test B_dual_C == expected_result
    end

    # C::Real
    C = test_value_2
    B_dual_C = dual(B, C)
    expected_result = test_value_1
    @test B_dual_C isa Scalar
    @test B_dual_C == expected_result
end
