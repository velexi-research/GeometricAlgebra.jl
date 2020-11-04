"""
Unit tests for the +(x, y) function

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

@testset "+(B::Pseudoscalar, C::Pseudoscalar)" begin
    test_dim = 10

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_plus_C = B + C
    @test B_plus_C isa Pseudoscalar
    @test B_plus_C == Pseudoscalar(test_dim, test_value_1 + test_value_2)
end

# ------ B::Scalar

@testset "+(B::Scalar, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == test_value_1 + test_value_2
end

@testset "+(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::Scalar, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Zero()

    @test B + C == B
end

@testset "+(B::Scalar, C::Real)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == test_value_1 + test_value_2
end

# ------ B::One

@testset "+(B::One, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::One, C::One)" begin
    B = One()
    C = One()

    C_plus_B = C + B
    @test B + C isa Scalar
    @test B + C == 2
end

@testset "+(B::One, C::Zero)" begin
    B = One()
    C = Zero()
    @test isone(B + C)
end

@testset "+(B::One, C::Real)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == 1 + test_value
end

# ------ B::Zero

@testset "+(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    @test B + C == C
end

@testset "+(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    @test isone(B + C)
end

@testset "+(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()
    @test iszero(B + C)
end

@testset "+(B::Zero, C::Real)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == C
end

# ------ B::Real

@testset "+(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == test_value_1 + test_value_2
end

@testset "+(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::Real, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    B_plus_C = B + C
    @test B_plus_C isa Scalar
    @test B_plus_C == B
end
