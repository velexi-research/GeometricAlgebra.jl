"""
Unit tests for the -(x, y) function

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

# ------ B::Pseudoscalar

@testset "-(B::Pseudoscalar, C::Pseudoscalar)" begin
    test_dim = 10

    # B != C
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_minus_C = B - C
    expected_result = Pseudoscalar(test_dim, test_value_1 - test_value_2)
    @test B_minus_C isa Pseudoscalar
    @test B_minus_C == expected_result

    # B == C
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_1)
    @test iszero(B - C)
end

# ------ B::Scalar

@testset "-(B::Scalar, C::Scalar)" begin
    # B != C
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = Scalar(test_value_1)
    C = Scalar(test_value_1)
    @test iszero(B - C)
end

@testset "-(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    B_minus_C = B - C
    expected_result = test_value - 1
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result
end

@testset "-(B::Scalar, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Zero()

    @test B - C === B
end

@testset "-(B::Scalar, C::Real)" begin
    # B != C
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = Scalar(test_value_1)
    C = test_value_1
    @test iszero(B - C)
end

# ------ B::One

@testset "-(B::One, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result
end

@testset "-(B::One, C::One)" begin
    B = One()
    C = One()
    @test iszero(B - C)
end

@testset "-(B::One, C::Zero)" begin
    B = One()
    C = Zero()
    @test isone(B - C)
end

@testset "-(B::One, C::Real)" begin
    # B != C
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = One()
    C = 1
    @test iszero(B - C)
end

# ------ B::Zero

@testset "-(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value
end

@testset "-(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -1
end

@testset "-(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()
    @test iszero(B - C)
end

@testset "-(B::Zero, C::Real)" begin
    # B != C
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == -test_value

    # B == C
    B = Zero()
    C = 0
    @test iszero(B - C)
end

# ------ B::Real

@testset "-(B::Real, C::Scalar)" begin
    # B != C
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_minus_C = B - C
    expected_result = test_value_1 - test_value_2
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = test_value_1
    C = Scalar(test_value_1)
    @test iszero(B - C)
end

@testset "-(B::Real, C::One)" begin
    # B != C
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_minus_C = B - C
    expected_result = test_value - 1
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    # B == C
    B = 1
    C = One()
    @test iszero(B - C)
end

@testset "-(B::Real, C::Zero)" begin
    # B != C
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C == test_value

    # B == C
    B = 0
    C = Zero()
    @test iszero(B - C)
end
