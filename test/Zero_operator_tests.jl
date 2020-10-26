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

@testset "-(B)" begin  # from AbstractMultivector interface
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        @test -B === B
    end
end

@testset "dual(B)" begin  # from AbstractMultivector interface
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        expected_message = "The dual of Zero is not well-defined"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "reciprocal(B)" begin  # from AbstractBlade interface
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar{precision_type}(Inf)
    end
end

@testset "+(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    @test B + C === Zero()

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    expected_result = One()
    @test B + C === expected_result
    @test C + B === expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value)

    @test B + C === C
    @test C + B === C

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

@testset "-(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    @test B - C === Zero()

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -1

    @test C - B === One()

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value

    @test C - B === C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value

    C_minus_B = C - B
    @test C_minus_B isa Scalar
    @test C_minus_B == test_value
end

@testset "*(B, C)" begin
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

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value)

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    expected_result = Zero()
    @test B * C === expected_result
    @test C * B === expected_result
end

@testset "/(B, C)" begin
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

    # B::Zero, C::One
    # B::One, C::Zero
    B = Zero()
    C = One()

    B_slash_C = B / C
    expected_result = Zero()
    @test B_slash_C === expected_result

    C_slash_B = C / B
    expected_result = Inf
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::Zero, C::Scalar
    # B::Scalar, C::Zero
    B = Zero()
    C = Scalar(test_value)

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
