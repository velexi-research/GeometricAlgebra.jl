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

# --- File inclusions

# Test utilities
include("test_utils.jl")

# --- Tests

# ------ M::Multivector

@testset "-(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "-(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "-(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::Scalar)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::One)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::Zero)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))
    C = Zero()
    @test B - C === B
end

@testset "-(B::Blade, C::Real)" begin
    @test_skip 1
end

@testset "-(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "-(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::Pseudoscalar, C::Pseudoscalar)" begin
    test_dim = 10

    # B != C
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_minus_C = B - C
    @test B_minus_C isa Pseudoscalar
    @test B_minus_C.dim == test_dim
    @test B_minus_C.value == test_value_1 - test_value_2

    # B == C
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_1)
    @test iszero(B - C)
end

@testset "-(B::Pseudoscalar, C::Scalar)" begin
    @test_skip 1
end

@testset "-(B::Pseudoscalar, C::One)" begin
    @test_skip 1
end

@testset "-(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0

    B = Pseudoscalar(test_dim, test_value)
    C = Zero()
    @test B - C === B
end

@testset "-(B::Pseudoscalar, C::Real)" begin
    @test_skip 1
end

@testset "-(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "-(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Scalar, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::Scalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(B::Scalar, C::Scalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    # --- Tests

    # B != C, C != B - 1
    test_value_2 = abs(test_value_1) + rand() + 2
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value_1 - test_value_2

    # B == C
    C = Scalar(test_value_1)
    @test iszero(B - C)

    # C == B + 1
    C = Scalar(test_value_1 - 1)
    
    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Scalar, C::One)" begin
    # --- Preparations

    C = One()

    # --- Tests

    # B != 2
    test_value = get_random_value(3)  # add 3 to keep value away from 0, 1, and 2
    B = Scalar(test_value)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value - 1

    # B == 2
    B = Scalar(2)

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Scalar, C::Zero)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = Zero()

    @test B - C === B
end

@testset "-(B::Scalar, C::Real)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    # --- Tests

    # B != C, C != B - 1
    test_value_2 = abs(test_value_1) + rand() + 2
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value_1 - test_value_2

    # B == C
    C = test_value_1
    @test iszero(B - C)

    # C == B - 1
    C= test_value_1 - 1

    B_minus_C = B - C
    @test B_minus_C isa One    
end

@testset "-(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "-(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::One, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(B::One, C::Scalar)" begin
    B = One()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == 1 - test_value
end

@testset "-(B::One, C::One)" begin
    B = One()
    C = One()
    @test iszero(B - C)
end

@testset "-(B::One, C::Zero)" begin
    B = One()
    C = Zero()

    B_minus_C = B - C
    @test B_minus_C === B
end

@testset "-(B::One, C::Real)" begin
    # B != C, C != B - 1
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == 1 - test_value

    # B == C
    B = One()
    C = 1
    @test iszero(B - C)

    # C == B - 1
    C = 0

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "-(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Zero, C::Blade)" begin
    test_dim = 10
    B = Zero()
    C = Blade(rand(test_dim, 3))
    @test B - C == -C
end

@testset "-(B::Zero, C::Pseudoscalar)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0

    B = Zero()
    C = Pseudoscalar(test_dim, test_value)
    @test B - C == -C
end

@testset "-(B::Zero, C::Scalar)" begin
    # --- Preparations

    B = Zero()

    # --- Tests

    # C != -1
    test_value = get_random_value(2)  # add 2 to keep value away from 0, 1, and -1
    C = Scalar(test_value)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == -test_value

    # C == -1
    C = Scalar(-1)

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == -1
end

@testset "-(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()
    
    B_minus_C = B - C
    @test B_minus_C === B
end

@testset "-(B::Zero, C::Real)" begin
    # --- Preparations

    B = Zero()

    # --- Tests

    # B != C, C != -1
    test_value = get_random_value(2)  # add 2 to keep value away from 0, 1, and -1
    C = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == -test_value

    # B == C
    B = Zero()
    C = 0
    @test iszero(B - C)

    # C == -1
    C = -1

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Zero, C::Vector)" begin
    B = Zero()
    C = rand(5)
    @test B - C == -C
end

# ------ B::Real

@testset "-(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Real, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::Real, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(B::Real, C::Scalar)" begin
    # --- Preparation

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value_1

    # --- Tests

    # B != C, C != B - 1
    test_value_2 = abs(test_value_1) + rand() + 2
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value_1 - test_value_2

    # B == C
    C = Scalar(test_value_1)
    @test iszero(B - C)

    # C == B - 1
    C = Scalar(test_value_1 - 1)

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Real, C::One)" begin
    # --- Preparations

    C = One()

    # --- Tests

    # B != C, B != C + 1
    test_value = get_random_value(3)  # add 2 to keep value away from 0, 1, and 2
    B = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value - 1

    # B == C
    B = 1
    @test iszero(B - C)

    # B = C + 1
    B = 2

    B_minus_C = B - C
    @test B_minus_C isa One
end

@testset "-(B::Real, C::Zero)" begin
    # --- Preparations

    C = Zero()

    # --- Tests

    # B != C, B != C + 1
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value

    B_minus_C = B - C
    @test B_minus_C isa Scalar
    @test B_minus_C.value == test_value

    # B == C
    B = 0
    C = Zero()
    @test iszero(B - C)

    # B == C + 1
    B = 1
    
    B_minus_C = B - C
    @test B_minus_C isa One
end

# ------ B::Vector

@testset "-(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "-(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "-(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "-(B::Vector, C::Scalar)" begin
    @test_skip 1
end

@testset "-(B::Vector, C::One)" begin
    @test_skip 1
end

@testset "-(B::Vector, C::Zero)" begin
    B = rand(5)
    C = Zero()
    @test B - C === B
end
