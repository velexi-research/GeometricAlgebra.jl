#   Copyright (c) 2020-2022 Velexi Corporation
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
Unit tests for the *(x, y) function
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

@testset "*(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "*(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "*(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "*(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "*(B::Blade, C::Scalar)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    @test B * C ≈ Blade(basis(B), volume=value(C) * volume(B))
end

@testset "*(B::Blade, C::One)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))
    C = One()

    @test B * C === B
end

@testset "*(B::Blade, C::Zero)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))
    C = Zero()

    @test B * C === C
end

@testset "*(B::Blade, C::Real)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = test_value

    @test B * C ≈ Blade(basis(B), volume=C * volume(B))
end

@testset "*(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "*(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

@testset "*(B::Pseudoscalar, C::Pseudoscalar)" begin
    # B != 1 / C
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    test_value_2 = abs(1 / test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_times_C = B * C
        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_times_C isa Scalar
        @test B_times_C.value == expected_result
    end

    # B == 1 / C
    # hard-coded a value for test_value_1 to avoid floating-point arithmetic errors
    # when computing the inverse of the value, which could result in the product
    # being different from one
    test_value_1 = 2

    for test_dim in 5:8
        test_value_2 = mod(test_dim, 4) < 2 ?
            1 / test_value_1 :
            -1 / test_value_1
        
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_times_C = B * C
        @test B_times_C isa One
    end

    # B ≈ 1 / C
    test_value_1 = 0.7621865484887302

    for test_dim in 5:8
        test_value_2 = mod(test_dim, 4) < 2 ?
            1 / test_value_1 :
            -1 / test_value_1
        
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_times_C = B * C
        @test B_times_C ≈ 1
    end
end

@testset "*(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C isa Pseudoscalar
    @test B_times_C.dim == test_dim
    @test B_times_C.value == test_value_1 * test_value_2
end

@testset "*(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value)

    C = One()

    @test B * C === B
end

@testset "*(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value)

    C = Zero()

    @test B * C === C
end

@testset "*(B::Pseudoscalar, C::Real)" begin
    # --- Preparation

    test_dim = 10
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = Pseudoscalar(test_dim, test_value_1)

    # --- Tests

    # C != 0
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0
    C = test_value_2

    B_times_C = B * C
    @test B_times_C isa Pseudoscalar
    @test B_times_C.dim == test_dim
    @test B_times_C.value == test_value_1 * test_value_2

    # C == 0
    C = 0

    B_times_C = B * C
    @test B_times_C isa Zero
end

@testset "*(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "*(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Scalar, C::Blade)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    @test B * C ≈ Blade(basis(C), volume=value(B) * volume(C))
end

@testset "*(B::Scalar, C::Pseudoscalar)" begin
    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    test_dim = 10
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0
    C = Pseudoscalar(test_dim, test_value_2)

    B_times_C = B * C
    @test B_times_C isa Pseudoscalar
    @test B_times_C.dim == test_dim
    @test B_times_C.value == test_value_1 * test_value_2
end

@testset "*(B::Scalar, C::Scalar)" begin
    # C != 1 / B
    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    test_value_2 = abs(1 / test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C isa Scalar
    @test B_times_C.value == test_value_1 * test_value_2

    # C == 1 / B
    # hard-coded a value for test_value_1 to avoid floating-point arithmetic errors
    # when computing the inverse of the value, which could result in the product
    # being different from one
    test_value_1 = 2
    B = Scalar(test_value_1)

    test_value_2 = 1 / test_value_1
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C isa One

    # C ≈ 1 / B
    test_value_1 = 0.7621865484887302
    B = Scalar(test_value_1)

    test_value_2 = 1 / test_value_1
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C ≈ 1
end

@testset "*(B::Scalar, C::One)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = One()

    @test B * C === B
end

@testset "*(B::Scalar, C::Zero)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = Zero()

    @test B * C === C
end

@testset "*(B::Scalar, C::Real)" begin
    # C != 0, C != 1 / B
    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    test_value_2 = abs(1 / test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_times_C = B * C
    @test B_times_C isa Scalar
    @test B_times_C.value == test_value_1 * test_value_2

    # C == 0
    C = 0

    B_times_C = B * C
    @test B_times_C isa Zero

    # C == 1 / B
    # hard-coded a value for test_value_1 to avoid floating-point arithmetic errors
    # when computing the inverse of the value, which could result in the product
    # being different from one
    test_value_1 = 2
    B = Scalar(test_value_1)

    test_value_2 = 1 / test_value_1
    C = test_value_2

    B_times_C = B * C
    @test B_times_C isa One

    # C ≈ 1 / B
    test_value_1 = 0.7621865484887302
    B = Scalar(test_value_1)

    test_value_2 = 1 / test_value_1
    C = test_value_2

    B_times_C = B * C
    @test B_times_C ≈ 1
end

@testset "*(B::Scalar, C::Vector)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    # C is not zero
    C = rand(5)

    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C ≈ test_value * Blade(C)

    # C is zero
    C = [0; 0; 0; 0; 0]
    @test iszero(B * C)
end

# ------ B::One

@testset "*(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::One, C::Blade)" begin
    B = One()
    test_dim = 10
    C = Blade(rand(test_dim, 3))

    @test B * C === C
end

@testset "*(B::One, C::Pseudoscalar)" begin
    B = One()

    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    C = Pseudoscalar(test_dim, test_value)

    @test B * C === C
end

@testset "*(B::One, C::Scalar)" begin
    B = One()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    @test B * C === C
end

@testset "*(B::One, C::One)" begin
    B = One()
    C = One()
    
    B_times_C = B * C
    @test B_times_C === B
end

@testset "*(B::One, C::Zero)" begin
    B = One()
    C = Zero()

    B_times_C = B * C
    @test B_times_C === C
end

@testset "*(B::One, C::Real)" begin
    # --- Preparations
    
    B = One()

    # --- Tests

    # C != 0, C != 1
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = test_value

    B_times_C = B * C
    @test B_times_C isa Scalar
    @test B_times_C.value == test_value

    # C == 0
    C = 0

    B_times_C = B * C
    @test B_times_C isa Zero

    # C == 1
    C = 1

    B_times_C = B * C
    @test B_times_C isa One
end

@testset "*(B::One, C::Vector)" begin
    C = One()

    # B is not zero
    B = rand(5)
    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C == Blade(B)

    # B is zero
    B = [0; 0; 0; 0; 0]
    @test iszero(B * C)
end

# ------ B::Zero

@testset "*(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Zero, C::Blade)" begin
    B = Zero()
    test_dim = 10
    C = Blade(rand(test_dim, 3))

    @test B * C === B
end

@testset "*(B::Zero, C::Pseudoscalar)" begin
    B = Zero()

    test_dim = 10
    test_value = get_random_value(1)  # add 1 to keep value away from 0
    C = Pseudoscalar(test_dim, test_value)

    @test B * C === B
end

@testset "*(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    @test B * C === B
end

@testset "*(B::Zero, C::One)" begin
    B = Zero()
    C = One()

    B_times_C = B * C
    @test B_times_C === B
end

@testset "*(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()

    B_times_C = B * C
    @test B_times_C === B
end

@testset "*(B::Zero, C::Real)" begin
    B = Zero()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = test_value

    @test B * C === B
end

@testset "*(B::Zero, C::Vector)" begin
    B = Zero()

    # C is not zero
    C = rand(5)
    @test B * C === B

    # C is zero
    C = [0; 0; 0; 0; 0]
    @test B * C === B
end

# ------ B::Real

@testset "*(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Real, C::Blade)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    @test B * C ≈ Blade(basis(C), volume=B * volume(C))
end

@testset "*(B::Real, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10
    test_value_2 = get_random_value(1)  # add 1 to keep value away from 0
    C = Pseudoscalar(test_dim, test_value_2)

    # --- Tests

    # B != 0
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = test_value_1

    B_times_C = B * C
    @test B_times_C isa Pseudoscalar
    @test B_times_C.dim == test_dim
    @test B_times_C.value == test_value_1 * test_value_2

    # B == 0
    B = 0

    B_times_C = B * C
    @test B_times_C isa Zero
end

@testset "*(B::Real, C::Scalar)" begin
    # B != 0, B != 1 / C
    test_value_1 = get_random_value(1)  # add 1 to keep value away from 0
    B = test_value_1

    test_value_2 = abs(1 / test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C isa Scalar
    @test B_times_C.value == test_value_1 * test_value_2

    # B == 0
    B = 0

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    B_times_C = B * C
    @test B_times_C isa Zero

    # B == 1 / C
    # hard-coded a value for test_value_1 to avoid floating-point arithmetic errors
    # when computing the inverse of the value, which could result in the product
    # being different from one
    test_value_1 = 2
    B = test_value_1

    test_value_2 = 1 / test_value_1
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C isa One

    # C ≈ 1 / B
    test_value_1 = 0.7621865484887302
    B = test_value_1

    test_value_2 = 1 / test_value_1
    C = Scalar(test_value_2)

    B_times_C = B * C
    @test B_times_C ≈ 1
end

@testset "*(B::Real, C::One)" begin
    # --- Preparations

    C = One()

    # --- Tests

    # B != 0, B != 1
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value

    B_times_C = B * C
    @test B_times_C isa Scalar
    @test B_times_C.value == test_value

    # B == 0
    B = 0

    B_times_C = B * C
    @test B_times_C isa Zero

    # B == 1
    B = 1

    B_times_C = B * C
    @test B_times_C isa One
end

@testset "*(B::Real, C::Zero)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = test_value

    C = Zero()

    @test iszero(B * C)
end

# ------ B::Vector

@testset "*(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "*(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "*(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "*(B::Vector, C::Scalar)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    # B is not zero
    B = rand(5)

    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C ≈ test_value * Blade(B)

    # B is zero
    B = [0; 0; 0; 0; 0]
    @test iszero(B * C)
end

@testset "*(B::Vector, C::One)" begin
    B = One()

    # C is not zero
    C = rand(5)
    B_times_C = B * C
    @test B_times_C isa Blade
    @test B_times_C == Blade(C)

    # C is zero
    C = [0; 0; 0; 0; 0]
    @test iszero(B * C)
end

@testset "*(B::Vector, C::Zero)" begin
    C = Zero()

    # B is not zero
    B = rand(5)
    @test B * C === C

    # B is zero
    B = [0; 0; 0; 0; 0]
    @test B * C === C
end
