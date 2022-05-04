#   Copyright 2020 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Unit tests for the /(x, y) function
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

@testset "/(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "/(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "/(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::Scalar)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::One)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))
    C = One()

    @test B / C === B
end

@testset "/(B::Blade, C::Zero)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::Real)" begin
    @test_skip 1
end

@testset "/(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "/(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Pseudoscalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10

    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    # --- Tests
    
    # B != C, C.value != Inf
    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == test_value_1 / test_value_2

    # B == C, C.value != Inf
    C = Pseudoscalar(test_dim, test_value_1)

    B_slash_C = B / C
    @test B_slash_C isa One

    # B != C, C.value == Inf
    test_value_2 = Inf
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == C, C.value == Inf
    test_value = Inf
    B = Pseudoscalar(test_dim, test_value)
    C = Pseudoscalar(test_dim, test_value)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end

@testset "/(B::Pseudoscalar, C::Scalar)" begin
    # --- Preparations

    test_dim = 10
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    # --- Tests

    # B != Inf, C != Inf
    test_value_2 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test B_slash_C.value == test_value_1 / test_value_2

    # B != Inf, C == Inf
    test_value_2 = Inf
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero
    
    # B == Inf, C == Inf
    test_value = Inf
    B = Pseudoscalar(test_dim, test_value)
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test isnan(B_slash_C.value)
end

@testset "/(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value)

    C = One()

    @test B / C === B
end

@testset "/(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value)

    C = Zero()

    expected_result = sign(B) > 0 ? Inf : -Inf

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test B_slash_C.value == expected_result
end

@testset "/(B::Pseudoscalar, C::Real)" begin
    # --- Preparations

    test_dim = 10
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    # --- Tests

    # B.value != Inf, C != Inf
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0
    C = test_value_2

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test B_slash_C.value == test_value_1 / test_value_2

    # B.value != Inf, C == Inf
    C = Inf

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B.value == Inf, C == Inf
    B = Pseudoscalar(test_dim, Inf)
    C = Inf

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test isnan(B_slash_C.value)
end

@testset "/(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "/(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Scalar, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Scalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    # --- Tests

    # B != Inf, C.value != Inf
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0

    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_2)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_1 / test_value_2) :
            Pseudoscalar(test_dim, -test_value_1 / test_value_2)

        B_slash_C = B / C
        @test B_slash_C isa Pseudoscalar
        @test B_slash_C.dim == expected_result.dim
        @test B_slash_C.value == expected_result.value
    end

    # B != Inf, C.value == Inf
    test_dim = 10
    test_value_2 = Inf
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == Inf, C.value == Inf
    test_dim = 10
    test_value = Inf
    B = Scalar(test_value)
    C = Pseudoscalar(test_dim, test_value)

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test isnan(B_slash_C.value)
end

@testset "/(B::Scalar, C::Scalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    # --- Tests

    # B != C, B != Inf, C != Inf
    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == test_value_1 / test_value_2

    # B == C, B != Inf
    C = Scalar(test_value_1)

    B_slash_C = B / C
    @test B_slash_C isa One

    # B != Inf, C == Inf
    C = Scalar(Inf)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == Inf, C == Inf
    test_value = Inf
    B = Scalar(test_value)
    C = Scalar(test_value)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end

@testset "/(B::Scalar, C::One)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = One()

    B_slash_C = B / C
    @test B_slash_C === B
end

@testset "/(B::Scalar, C::Zero)" begin
    # --- Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1

    C = Zero()

    # --- Tests

    # B > 0
    B = Scalar(abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == Inf

    # B < 0
    B = Scalar(-abs(test_value))
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == -Inf
end

@testset "/(B::Scalar, C::Real)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    # --- Tests

    # B != C, B != Inf, C != Inf, C != 0
    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == test_value_1 / test_value_2

    # B == C, B != Inf, C != 0
    C = test_value_1

    B_slash_C = B / C
    @test B_slash_C isa One

    # B != Inf, C == Inf
    C = Inf

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B > 0, B != Inf, C == 0
    B = Scalar(abs(test_value_1))
    C = 0

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == Inf

    # B < 0, B != -Inf, C == 0
    B = Scalar(-abs(test_value_1))
    C = 0

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == -Inf

    # B == Inf, C == Inf
    test_value = Inf
    B = Scalar(test_value)
    C = test_value

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end

@testset "/(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "/(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::One, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0

    B = One()

    # --- Tests

    # C.value != Inf
    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_1)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, 1/ test_value_1) :
            Pseudoscalar(test_dim, -1/ test_value_1)

        B_slash_C = B / C
        @test B_slash_C isa Pseudoscalar
        @test B_slash_C.dim == expected_result.dim
        @test B_slash_C.value == expected_result.value
    end

    # C.value == Inf
    test_dim = 10
    test_value_2 = Inf
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero
end

@testset "/(B::One, C::Scalar)" begin
    # C != Inf
    B = One()

    test_value = get_random_value(2)  # add 1 to keep value away from 0 and 1
    C = Scalar(test_value)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == 1 / test_value

    # C == Inf
    C = Scalar(Inf)

    B_slash_C = B / C
    @test B_slash_C isa Zero
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

    # B != C, C != Inf, C != 0
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = test_value

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == 1 / test_value

    # B == C
    C = 1

    B_slash_C = B / C
    @test B_slash_C isa One

    # C == Inf
    C = Inf

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # C == 0
    C = 0

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == Inf
end

@testset "/(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "/(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Zero, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Zero, C::Pseudoscalar)" begin
    B = Zero()

    test_dim = 10
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    C = Pseudoscalar(test_dim, test_value_1)

    @test iszero(B / C)
end

@testset "/(B::Zero, C::Scalar)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
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
    @test isnan(B_slash_C.value)
end

@testset "/(B::Zero, C::Real)" begin
    # --- Preparations

    B = Zero()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1

    # --- Tests

    # C != 0
    C = test_value

    @test iszero(B / C)

    # C == 0
    C = 0

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end

@testset "/(B::Zero, C::Vector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "/(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Real, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Real, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value_1

    # --- Tests

    # B != Inf, B != 0, C.value != Inf
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0

    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value_2)

        expected_result = mod(test_dim, 4) < 2 ?
            Pseudoscalar(test_dim, test_value_1 / test_value_2) :
            Pseudoscalar(test_dim, -test_value_1 / test_value_2)

        B_slash_C = B / C
        @test B_slash_C isa Pseudoscalar
        @test B_slash_C.dim == expected_result.dim
        @test B_slash_C.value == expected_result.value
    end

    # B != Inf, B != 0, C.value == Inf
    test_dim = 10
    test_value_2 = Inf
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == 0
    B = 0

    test_dim = 10
    C = Pseudoscalar(test_dim, test_value_1)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == Inf, C.value == Inf
    test_value = Inf
    test_dim = 10
    B = test_value
    C = Pseudoscalar(test_dim, test_value)

    B_slash_C = B / C
    @test B_slash_C isa Pseudoscalar
    @test B_slash_C.dim == test_dim
    @test isnan(B_slash_C.value)
end

@testset "/(B::Real, C::Scalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value_1

    # --- Tests

    # B != C, B != Inf, C != Inf
    test_value_2 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == test_value_1 / test_value_2

    # B == C, B != Inf
    C = Scalar(test_value_1)

    B_slash_C = B / C
    @test B_slash_C isa One

    # B != Inf, C == Inf
    test_value_2 = Inf
    C = Scalar(test_value_2)

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == 0
    B = 0

    B_slash_C = B / C
    @test B_slash_C isa Zero

    # B == Inf, C == Inf
    test_value = Inf
    B = test_value
    C = Scalar(test_value)

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end

@testset "/(B::Real, C::One)" begin
    # --- Preparations

    C = One()

    # --- Tests
    
    # B != C, B != 0
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value

    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C.value == test_value

    # B == C
    B = 1

    B_slash_C = B / C
    @test B_slash_C isa One

    # B == 0
    B = 0

    B_slash_C = B / C
    @test B_slash_C isa Zero
end

@testset "/(B::Real, C::Zero)" begin
    # --- Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1

    C = Zero()

    # --- Tests

    # B > 0
    B = abs(test_value)
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == Inf

    # B = 0
    B = 0
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)

    # B < 0
    B = -abs(test_value)
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test B_slash_C == -Inf
end

# ------ B::Vector

@testset "/(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "/(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "/(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "/(B::Vector, C::Scalar)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    # B is not zero
    B = rand(5)

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B) / test_value

    # B is zero
    B = [0; 0; 0; 0; 0]
    @test iszero(B / C);
end

@testset "/(B::Vector, C::One)" begin
    C = One()

    # B is not zero
    B = rand(5)

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B)

    # B is zero
    B = [0; 0; 0; 0; 0]
    @test iszero(B / C)
end

@testset "/(B::Vector, C::Zero)" begin
    C = Zero()

    # B is not zero
    B = rand(5)

    B_slash_C = B / C
    @test B_slash_C isa Blade
    @test B_slash_C == Blade(B) / 0

    # B is zero
    B = [0; 0; 0; 0; 0]
    
    B_slash_C = B / C
    @test B_slash_C isa Scalar
    @test isnan(B_slash_C.value)
end
