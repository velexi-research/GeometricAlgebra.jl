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

@testset "-(B)" begin  # from AbstractMultivector interface
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        minus_B = -B
        @test minus_B isa Scalar{precision_type}
        @test minus_B == Scalar{precision_type}(-1)
    end
end

@testset "dual(B)" begin  # from AbstractMultivector interface
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)
            @test dual_B isa Pseudoscalar{precision_type}

            expected_dual = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1) :
                Pseudoscalar{precision_type}(test_dim, -1)
            @test dual_B == expected_dual
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "reciprocal(B)" begin  # from AbstractBlade interface
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B === B
    end
end

@testset "+(B, C)" begin
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

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_plus_C = B + C
    expected_result = 1 + test_value
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

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

@testset "-(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    expected_result = Zero()
    @test B - C === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
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
end

@testset "*(B, C)" begin
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

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    expected_result = C
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
end

@testset "/(B, C)" begin
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

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_slash_C = B / C
    expected_result = 1 / test_value
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

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
end
