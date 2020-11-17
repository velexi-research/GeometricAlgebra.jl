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

# ------ M::Multivector

@testset "+(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "+(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "+(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Scalar)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::One)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Zero)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Real)" begin
    @test_skip 1
end

@testset "+(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "+(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

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

@testset "+(B::Pseudoscalar, C::Scalar)" begin
    @test_skip 1
end

@testset "+(B::Pseudoscalar, C::One)" begin
    @test_skip 1
end

@testset "+(B::Pseudoscalar, C::Zero)" begin
    @test_skip 1
end

@testset "+(B::Pseudoscalar, C::Real)" begin
    @test_skip 1
end

@testset "+(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "+(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Scalar, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::Scalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(B::Scalar, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == test_value_1 + test_value_2
end

@testset "+(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
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
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == test_value_1 + test_value_2
end

@testset "+(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "+(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::One, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(B::One, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::One, C::One)" begin
    B = One()
    C = One()

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == 2
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
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "+(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Zero, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::Zero, C::Pseudoscalar)" begin
    @test_skip 1
end

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
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == C
end

@testset "+(B::Zero, C::Vector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "+(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Real, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::Real, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == test_value_1 + test_value_2
end

@testset "+(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == 1 + test_value
end

@testset "+(B::Real, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    B_plus_C = B + C
    @test B_plus_C isa AbstractScalar
    @test B_plus_C == B
end

# ------ B::Vector

@testset "+(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "+(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "+(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "+(B::Vector, C::Scalar)" begin
    @test_skip 1
end

@testset "+(B::Vector, C::One)" begin
    @test_skip 1
end

@testset "+(B::Vector, C::Zero)" begin
    @test_skip 1
end
