"""
Unit tests for methods defined for the Scalar type.

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

# ------ Unary operations

@testset "Scalar: -(B)" begin
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

@testset "Scalar: reverse(B)" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    # value != 0
    B = Scalar(test_value)
    reverse_B = reverse(B)
    @test reverse_B === B
    @test B * reverse_B ≈ norm(B)^2

    # value = Inf
    B = Scalar(Inf)
    @test reverse(B) === B

    # value = -Inf
    B = Scalar(-Inf)
    @test reverse(B) === B
end

@testset "Scalar: dual(B)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)

            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(test_value)) :
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(-test_value))
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_result
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "Scalar: reciprocal(B)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)

        # value > 0
        B = Scalar(abs(converted_value))
        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar(1 / abs(converted_value))
        @test B * reciprocal_B ≈ 1

        # value < 0
        negative_value = -(abs(converted_value))
        B = Scalar(negative_value)
        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar(1 / negative_value)
        @test B * reciprocal_B ≈ 1

        # value = Inf
        B = Scalar(precision_type(Inf))
        @test reciprocal(B) === Zero{precision_type}()

        # value = -Inf
        B = Scalar(precision_type(-Inf))
        @test reciprocal(B) === Zero{precision_type}()
    end
end

# ------ Binary operations

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

    # B::Vector, C::Scalar
    B = rand(5)
    C = Scalar(test_value_2)

    expected_result = Blade(B) / C

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == expected_result
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
    @test B ∧ C == B_wedge_C

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = B
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = B
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B

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
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B
end

@testset "Scalar: contractl(B, C), B < C" begin
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
end

@testset "Scalar: dual(B, C)" begin
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

@testset "Scalar: proj(B, C)" begin
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

    B_proj_C = proj(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    # --- B::Scalar, C::One
    #     B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_proj_C = proj(B, C)
    expected_result = B
    @test B_proj_C == expected_result

    C_proj_B = proj(C, B)
    expected_result = C
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Zero
    #     B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_proj_C = proj(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Real
    #     B::Real, C::Scalar
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

    # --- B::Scalar, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::Scalar

    # return_blade == true
    B = Scalar(test_value_1)
    C = test_vector
    @test proj(B, C) == B
    @test proj(C, B) === zero(B)

    # return_blade == false
    @test proj(C, B, return_blade=false) == 0
    @test proj(B, C, return_blade=false) == value(B)
end

# ------ Comparison operators

@testset "Scalar: ==(B, C)" begin
    # Preparations
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    float64_or_bigfloat = (Float64, BigFloat)

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(B, C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            if precision_type1 == precision_type2
                @test B == C
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B == C
            else
                @test B != C
            end

            # value(B) != value(C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(2 * test_value))
            @test B != C
        end
    end

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    for precision_type in subtypes(AbstractFloat)
        # B::Scalar, C::AbstractFloat
        # B::AbstractFloat, C::Scalar
        for value_type in subtypes(AbstractFloat)
            B = Scalar(precision_type(test_value))

            # ==(B, C)
            if precision_type == value_type
                @test B == value_type(test_value)
                @test value_type(test_value) == B
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test B == value_type(test_value)
                @test value_type(test_value) == B
            else
                @test B != value_type(test_value)
                @test value_type(test_value) != B
            end

            # !=(x,y)
            @test B != value_type(2 * test_value)
            @test value_type(2 * test_value) != B
        end

        # B::Scalar, C::Integer
        # B::Integer, C::Scalar
        int_value::Int = 3
        for value_type in subtypes(Signed)
            B = Scalar(precision_type(int_value))

            # ==(B, C)
            @test B == value_type(int_value)
            @test value_type(int_value) == B

            # !=(x,y)
            @test B != value_type(2 * int_value)
            @test value_type(2 * int_value) != B
        end

        for value_type in subtypes(Unsigned)
            B = Scalar(precision_type(int_value))

            # ==(B, C)
            @test B == value_type(int_value)
            @test value_type(int_value) == B

            # !=(x,y)
            @test B != value_type(2 * int_value)
            @test value_type(2 * int_value) != B
        end

        # Bool
        B = Scalar(precision_type(true))
        @test B == 1
        @test 1 == B
        @test B != 0
        @test 0 != B

        B = Scalar(precision_type(false))
        @test B == 0
        @test 0 == B
        @test B != 1
        @test 1 != B
    end
end

@testset "Scalar: ≈(B, C)" begin
    # Preparations
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            @test B ≈ C
        end
    end

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)
        B = Scalar(converted_value)
        @test B ≈ converted_value
        @test converted_value ≈ B
    end
end

# ------ Utility methods

@testset "Scalar: convert(S)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Exercise functionality and check results
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            converted_test_value = precision_type_src(test_value)
            S = Scalar{precision_type_src}(converted_test_value)
            S_converted = convert(AbstractScalar{precision_type_converted}, S)
            @test S_converted isa Scalar{precision_type_converted}
            if precision_type_src == precision_type_converted
                @test S_converted === S
            else
                @test S_converted !== S
                @test S_converted ≈ S
            end
        end
    end
end
