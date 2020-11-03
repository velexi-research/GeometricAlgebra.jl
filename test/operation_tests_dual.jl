"""
Unit tests for the dual(x, y) function

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
import LinearAlgebra
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

@testset "dual(B, C): B::Blade, C::Blade" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 20

    C = Blade(rand(test_dim, 7))

    # --- Exercise functionality and check results

    # ------ dim(B) == dim(C), grade(B) < grade(C)

    for grade_B in 2:5
        B_vectors = basis(C)[:, 1:grade_B] * rand(grade_B, grade_B)
        B = Blade(B_vectors)

        dual_B = dual(B, C)

        # --- Check dim, grade, and norm

        @test dim(dual_B) == test_dim
        @test grade(dual_B) == grade(C) - grade_B
        @test norm(dual_B) == norm(B)

        # --- Check that B and dual(B) are orthogonal complements

        @test LinearAlgebra.norm(transpose(basis(dual_B)) * basis(B)) <
            10 * eps(Float64)

        # --- Check sign(dual(B, C))

        # Compute sign of I_C formed from basis(B) and basis(dual(B, C))
        sign_Q = sign(LinearAlgebra.det(
            transpose(hcat(basis(B), basis(dual_B))) * basis(C)))

        # Compute expected_sign
        expected_sign = sign(B) * sign_Q

        # Account for sign of I_C^{-1} relative to I_C
        if mod(grade(C), 4) >= 2
            expected_sign = -expected_sign
        end

        # Account for reversals required to eliminate B
        if mod(grade(B), 4) >= 2
            expected_sign = -expected_sign
        end

        @test sign(dual_B) == expected_sign

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) ≈ B
        else
            @test dual(dual_B, C) ≈ -B
        end
    end

    # ------ dim(B) == dim(C), grade(B) == grade(C)

    coefficients = rand(grade(C), grade(C))
    B_vectors = basis(C) * coefficients
    B = Blade(B_vectors)

    relative_orientation =
        sign(LinearAlgebra.det(transpose(basis(B)) * basis(C)))
    expected_result = mod(grade(B), 4) < 2 ?
        Scalar(relative_orientation * volume(B)) :
        Scalar(-relative_orientation * volume(B))
    @test dual(B, C) ≈ expected_result

    # ------ Error cases

    for grade_B in 2:5
        # dim(B) != dim(C)
        B = Blade(rand(test_dim + 1, grade_B))
        @test_throws DimensionMismatch dual(B, C)

        # `B` not contained in `C`
        B = Blade(rand(test_dim, grade_B))
        while LinearAlgebra.norm(reject(basis(B), C)) < sqrt(eps(Float64))
            B = Blade(rand(test_dim, grade_B))
        end
        @test_throws ErrorException dual(B, C)
    end
end

@testset "dual(B, C): B or C isa Blade" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Blade, C::Pseudoscalar
    #     B::Pseudoscalar, C::Blade

    # dim(B) == dim(C)
    for grade_B in 2:5
        B = Blade(rand(test_dim, grade_B))
        C = Pseudoscalar(test_dim, test_value_1)

        dual_B = dual(B, C)
        expected_result = dual(B)
        @test dual_B == expected_result

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) ≈ B
        else
            @test dual(dual_B, C) ≈ -B
        end

        dual_B = dual(C, B)
        expected_result = zero(C)
        @test dual_B === expected_result
    end

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Pseudoscalar(test_dim + 1, test_value_1)
    @test_throws DimensionMismatch dual(B, C)
    @test_throws DimensionMismatch dual(C, B)

    # --- B::Blade, C::Scalar
    #     B::Scalar, C::Blade

    B = Scalar(test_value_1)

    for test_grade in 5:8
        C = Blade(randn(test_dim, test_grade))
        B_dual_C = dual(B, C)
        expected_result = mod(test_grade, 4) < 2 ?
            Blade(C, volume=test_value_1) :
            Blade(C, volume=-test_value_1)
        @test B_dual_C isa Blade
        @test B_dual_C == expected_result

        C_dual_B = dual(C, B)
        expected_result = zero(C)
        @test C_dual_B === expected_result
    end

    # --- B::Blade, C::One
    #     B::One, C::Blade

    B = One()

    for test_grade in 5:8
        C = Blade(randn(test_dim, test_grade))

        B_dual_C = dual(B, C)
        expected_result = mod(test_grade, 4) < 2 ?
            Blade(C, volume=1, copy_basis=false) :
            Blade(C, volume=-1, copy_basis=false)
        @test B_dual_C isa Blade
        @test B_dual_C == expected_result

        C_dual_B = dual(C, B)
        expected_result = zero(C)
        @test C_dual_B === expected_result
    end
end

@testset "dual(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)
        dual_B = dual(B, C)
        expected_result = test_value_1
        @test dual_B isa Scalar
        @test dual_B == expected_result

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) == B
        else
            @test dual(dual_B, C) == -B
        end
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dual(B, C)

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    C = Scalar(test_value_2)

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) === expected_result

        expected_result = mod(dim_B, 4) < 2 ?
            Pseudoscalar(dim_B, test_value_2) :
            Pseudoscalar(dim_B, -test_value_2)
        @test dual(C, B) == expected_result
    end

    # --- B::Pseudoscalar, C::One
    #     B::One, C::Pseudoscalar

    C = One()

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result

        expected_result = mod(dim_B, 4) < 2 ?
            Pseudoscalar(dim_B, 1) :
            Pseudoscalar(dim_B, -1)
        @test dual(C, B) == expected_result
    end

    # --- B::Pseudoscalar, C::Real
    #     B::Real, C::Pseudoscalar

    C = test_value_2

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result
    end
end

@testset "dual(B, C): B or C isa Scalar" begin
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

@testset "dual(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    B = One()

    # B::One, C::One
    C = One()
    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    # B::One, C::Scalar
    # C::Scalar, B::One
    C = Scalar(test_value)

    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    C_dual_B = dual(C, B)
    expected_result = C
    @test C_dual_B == expected_result

    # B::One, C::Real
    # C::Real, C::One
    C = test_value

    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    C_dual_B = dual(C, B)
    expected_result = Scalar(C)
    @test C_dual_B == expected_result
end

@testset "dual(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- dual(B::Zero, C)

    B = Zero()
    expected_error = "The dual of Zero is not well-defined"

    # C::One
    C = One()
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Scalar
    C = Scalar(test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Real
    C = test_value
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Blade
    C = Blade(randn(4, 3))
    @test_throws ErrorException(expected_error) dual(B, C)

    # C::Pseudoscalar
    C = Pseudoscalar(5, test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # --- dual(B, C::Zero)

    C = Zero()

    expected_error = "The dual of anything relative to Zero is not well-defined"

    # B::Zero
    B = Zero()
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::One
    B = One()
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Scalar
    B = Scalar(test_value)
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Real
    B = test_value
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Blade
    B = Blade(randn(4, 3))
    @test_throws ErrorException(expected_error) dual(B, C)

    # B::Pseudoscalar
    B = Pseudoscalar(5, test_value)
    @test_throws ErrorException(expected_error) dual(B, C)
end
