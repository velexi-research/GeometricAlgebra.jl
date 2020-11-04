"""
Unit tests for the /(x, y) function

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

@testset "/(B::Pseudoscalar, C::Pseudoscalar)" begin
    test_dim = 10

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value_1 / test_value_2
end

@testset "/(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    @test B / C == Pseudoscalar(test_dim, test_value_1 / test_value_2)
end

@testset "/(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = One()

    @test B / C === B
end

@testset "/(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = Zero()

    expected_result = sign(B) > 0 ?
        Pseudoscalar(B, value=Inf) :
        Pseudoscalar(B, value=-Inf)

    @test B / C == expected_result
end

@testset "/(B::Pseudoscalar, C::Real)" begin
    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    @test B / C == Pseudoscalar(test_dim, test_value_1 / test_value_2)
end

# ------ B::Scalar

@testset "/(B::Scalar, C::Pseudoscalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_2)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_1 / test_value_2) :
            Pseudoscalar(test_dim, -test_value_1 / test_value_2)

        @test B / C == expected_result
    end
end

@testset "/(B::Scalar, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value_1 / test_value_2
end

@testset "/(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value
end

@testset "/(B::Scalar, C::Zero)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    C = Zero()

    # --- Tests

    # B > 0
    B = Scalar(abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == Inf

    # B < 0
    B = Scalar(-abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == -Inf
end

@testset "/(B::Scalar, C::Real)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value_1 / test_value_2
end

# ------ B::One

@testset "/(B::One, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    B = One()

    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, 1/ test_value) :
            Pseudoscalar(test_dim, -1/ test_value)

        @test B / C == expected_result
    end
end

@testset "/(B::One, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == 1 / test_value
end

@testset "/(B::One, C::One)" begin
    B = One()
    C = One()
    @test isone(B / C)
end

@testset "/(B::One, C::Zero)" begin
    B = One()
    C = Zero()

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == Inf
end

@testset "/(B::One, C::Real)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == 1 / test_value
end

# ------ B::Zero

@testset "/(B::Zero, C::Pseudoscalar)" begin
    B = Zero()

    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    C = Pseudoscalar(test_dim, test_value_1)

    @test iszero(B / C)
end

@testset "/(B::Zero, C::Scalar)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Zero()

    @test iszero(C / B)
end

@testset "/(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    @test iszero(B / C)
end

@testset "/(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(value(B_slash_C))
end

@testset "/(B::Zero, C::Real)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    @test iszero(B / C)
end

# ------ B::Real

@testset "/(B::Real, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_2)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_1 / test_value_2) :
            Pseudoscalar(test_dim, -test_value_1 / test_value_2)

        @test B / C == expected_result
    end
end

@testset "/(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value_1 / test_value_2
end

@testset "/(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == test_value
end

@testset "/(B::Real, C::Zero)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    C = Zero()

    # --- Tests

    # B > 0
    B = Scalar(abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == Inf

    # B < 0
    B = Scalar(-abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == -Inf
end

# ------ B::Vector

@testset "/(B::Vector, C::Scalar)" begin
    B = rand(5)

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B) / test_value
end

@testset "/(B::Vector, C::One)" begin
    B = rand(5)
    C = One()

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B)
end

@testset "/(B::Vector, C::Zero)" begin
    B = rand(5)
    C = Zero()

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B) / 0
end
