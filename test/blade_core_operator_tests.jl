"""
Unit tests for AbstractBlade binary operators.

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
import LinearAlgebra
using Test

# GeometricAlgebra.jl
using GeometricAlgebra


# --- *(B, C): scalar multiplication

@testset "*(B, C): B or C isa {Scalar, Real}" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # Blade vectors
    vectors = randn(test_dim, 3)

    # --- x::Real, B::Scalar
    #     B::Scalar, x::Real

    # Preparations
    x = test_value_1
    B = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, x * value(B))
    @test x * B == expected_result
    @test B * x == expected_result

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B) * value(C))
    @test B * C == expected_result

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x * B ≈ expected_result
    @test B * x ≈ expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Blade(vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test B * C ≈ expected_result
    @test C * B ≈ expected_result

    # --- x::Real, B::Pseudoscalar
    #     B::Pseudoscalar, x::Real

    # Preparations
    x = test_value_1
    B = Pseudoscalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Pseudoscalar(test_dim, x * value(B))
    @test x * B == expected_result
    @test B * x == expected_result

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Pseudoscalar(test_dim, value(B) * value(C))
    @test B * C == expected_result
    @test C * B == expected_result
end

# --- ∧(B, C), outer(B, C)

@testset "∧(B, C): B or C isa {Scalar, Real}" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # Blade vectors
    vectors = randn(test_dim, 3)

    # --- x::Real, B::Scalar
    #     B::Scalar, x::Real

    # Preparations
    x = test_value_1
    B = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, x * value(B))
    @test x ∧ B == expected_result
    @test B ∧ x == expected_result
    @test outer(x, B) == expected_result
    @test outer(B, x) == expected_result

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B) * value(C))
    @test B ∧ C == expected_result
    @test outer(B, C) == expected_result

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x ∧ B ≈ expected_result
    @test B ∧ x ≈ expected_result
    @test outer(x, B) ≈ expected_result
    @test outer(B, x) ≈ expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Blade(vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test B ∧ C ≈ expected_result
    @test C ∧ B ≈ expected_result
    @test outer(B, C) ≈ expected_result
    @test outer(C, B) ≈ expected_result

    # --- x::Real, B::Pseudoscalar
    #     B::Pseudoscalar, x::Real

    # Preparations
    x = test_value_1
    B = Pseudoscalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Pseudoscalar(test_dim, x * value(B))
    @test x ∧ B == expected_result
    @test B ∧ x == expected_result
    @test outer(x, B) == expected_result
    @test outer(B, x) == expected_result

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Pseudoscalar(test_dim, value(B) * value(C))
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result
    @test outer(B, C) == expected_result
    @test outer(C, B) == expected_result
end

@testset "∧(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # Blade vectors
    vectors = randn(test_dim, 3)

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = zero(B)
    @test B ∧ C == expected_result
    @test outer(B, C) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch B ∧ C

    # --- B::Pseudoscalar, C::Blade
    #     B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(vectors)

    expected_result = zero(B)
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result
    @test outer(B, C) == expected_result
    @test outer(C, B) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(vectors)
    @test_throws DimensionMismatch B ∧ C
    @test_throws DimensionMismatch C ∧ B
end

@testset "∧(B, C): B, C::{Blade, Vector}" begin
    # --- B::Blade, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vectors = hcat([0; 0; 3; 0; 0],
                     [0; 0; 0; 4; 0])
    C = Blade(C_vectors)

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C ≈ Blade(hcat(B_vectors, C_vectors))
    @test B_wedge_C == outer(B, C)

    # --- B::Blade, C::Vector
    #     B::Vector, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vector = [0; 0; 3; 0; 0]
    C = Blade(C_vector)

    expected_result = Blade(hcat(B_vectors, C_vector))

    # Exercise functionality
    B_wedge_C = B ∧ C_vector
    @test B_wedge_C ≈ expected_result
    @test B_wedge_C == outer(B, C)

    C_wedge_B = C_vector ∧ B
    @test C_wedge_B ≈ (-1)^(grade(B)) * expected_result
    @test C_wedge_B == outer(C, B)

    # --- B::Vector, C::Vector

    # Preparations
    B_vector = Vector{Float64}([0, 2, 0, 0, 0])
    B = Blade(B_vector)

    C_vector = [0; 0; 3; 0; 0]
    C = Blade(C_vector)

    expected_result = Blade(hcat(B_vector, C_vector))

    # Exercise functionality
    B_wedge_C = B_vector ∧ C_vector
    @test B_wedge_C ≈ expected_result
    @test B_wedge_C == outer(B, C)
end

# --- project(B, C)

@testset "project(B, C): B or C isa Vector" begin
    @test_skip 1
end

@testset "project(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # Blade vectors
    vectors = randn(test_dim, 3)

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B))
    @test project(B, C) == expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Blade(vectors)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B))
    @test project(B, C) == expected_result

    expected_result = zero(C)
    @test project(C, B) == expected_result

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B))
    @test project(B, C) == expected_result

    expected_result = zero(C)
    @test project(C, B) == expected_result
end

@testset "project(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # Blade vectors
    vectors = randn(test_dim, 3)

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = B
    @test project(B, C) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch project(B, C)

    # --- B::Pseudoscalar, C::Blade
    #     B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(vectors)

    expected_result = zero(B)
    @test project(B, C) == expected_result

    expected_result = C
    @test project(C, B) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(vectors)
    @test_throws DimensionMismatch project(B, C)
    @test_throws DimensionMismatch project(C, B)
end

@testset "project(B, C): B, C::Blade" begin
    @test_skip 1
end

@testset "dual(B) tests" begin
    # --- Preparations

    # Test values
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = 20

    # --- B::Scalar

    B = Scalar(test_dim, test_value)

    dual_B = dual(B)

    expected_result = mod(dim(B), 4) < 2 ?
        Pseudoscalar(test_dim, value(B)) : Pseudoscalar(test_dim, -value(B))
    @test dual_B == expected_result

    # --- B::Pseudoscalar

    B = Pseudoscalar(test_dim, test_value)

    dual_B = dual(B)

    expected_result = Scalar(test_dim, value(B))
    @test dual_B == expected_result

    # Verify dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
    if mod(dim(B), 4) < 2
        @test_skip dual(dual_B) ≈ B
    else
        @test_skip dual(dual_B) ≈ -B
    end

    # --- B::Blade

    for dim_B in 10:13
        for grade_B in 2:5
            B = Blade(randn(dim_B, grade_B))

            dual_B = dual(B)

            # --- Verify dim, grade, and norm

            @test dim(dual_B) == dim_B
            @test grade(dual_B) == dim_B - grade_B
            @test norm(dual_B) == norm(B)

            # --- Verify that B and dual(B) are orthogonal complements

            @test LinearAlgebra.norm(transpose(basis(dual_B)) * basis(B)) <
                10 * eps(Float64)

            # --- Verify sign(dual(B))

            # Compute sign of I formed from basis(B) and basis(dual(B))
            sign_Q = sign(LinearAlgebra.det(hcat(basis(B), basis(dual_B))))

            # Compute expected_sign
            expected_sign = sign(B) * sign_Q

            # Account for sign of I^{-1} relative to I
            if mod(dim(B), 4) >= 2
                expected_sign = -expected_sign
            end

            # Account for reversals required to eliminate B
            if mod(grade(B), 4) >= 2
                expected_sign = -expected_sign
            end

            @test sign(dual_B) == expected_sign

            # --- Verify mathematical properties of dual

            # dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
            if mod(dim(B), 4) < 2
                @test dual(dual_B) ≈ B
            else
                @test dual(dual_B) ≈ -B
            end
        end
    end
end

# --- dual(B, C)

@testset "dual(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Dimension of embedding space
    test_dim = 10

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)
    C = Scalar(test_dim, test_value_2)

    # Exercise functionality and check results
    expected_result = Scalar(test_dim, value(B))
    @test dual(B, C) == expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)

    # Exercise functionality and check results
    for grade_C in 1:4
        C = Blade(randn(test_dim, grade_C))

        expected_result = mod(grade(C), 4) < 2 ?
            Blade(C, volume=value(B)) :
            Blade(C, volume=-value(B))
        @test dual(B, C) ≈ expected_result

        expected_result = zero(C)
        @test dual(C, B) == expected_result
    end

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_dim, test_value_1)

    # Exercise functionality and check results
    for dim_C in test_dim:test_dim + 3
        C = Pseudoscalar(dim_C, test_value_2)

        expected_result = mod(grade(C), 4) < 2 ?
            Pseudoscalar(dim_C, value(B)) :
            Pseudoscalar(dim_C, -value(B))
        @test dual(B, C) == expected_result

        expected_result = zero(C)
        @test dual(C, B) == expected_result
    end
end

@testset "dual(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    test_dim = 10

    # --- B::Pseudoscalar, C::Pseudoscalar

    # ------ dim(B) == dim(C)

    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    dual_B = dual(B, C)

    expected_result = Scalar(test_dim, value(B))
    @test dual_B == expected_result

    # Verify dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
    if mod(grade(C), 4) < 2
        @test dual(dual_B, C) ≈ B
    else
        @test dual(dual_B, C) ≈ -B
    end

    # ------ Error cases

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dual(B, C)

    # --- B::Blade, C::Pseudoscalar
    #     B::Pseudoscalar, C::Blade

    # ------ dim(B) == dim(C)

    B = Blade(randn(test_dim, 3))
    C = Pseudoscalar(test_dim, test_value_1)

    expected_result = zero(B)
    @test dual(C, B) == expected_result

    for grade_B in 2:5
        B = Blade(randn(test_dim, grade_B))
        C = Pseudoscalar(test_dim, test_value_1)

        dual_B = dual(B, C)
        @test dual_B == dual(B)

        # Verify dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) ≈ B
        else
            @test dual(dual_B, C) ≈ -B
        end
    end

    # ------ Error cases

    # dim(B) != dim(C)
    B = Blade(randn(test_dim, 3))
    C = Pseudoscalar(test_dim + 1, test_value_1)
    @test_throws DimensionMismatch dual(B, C)
    @test_throws DimensionMismatch dual(C, B)
end

@testset "dual(B, C): B, C::Blade" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 20

    C = Blade(randn(test_dim, 7))

    # --- Exercise functionality and check results

    # ------ dim(B) == dim(C), grade(B) < grade(C)

    for grade_B in 2:5
        B_vectors = basis(C)[:, 1:grade_B] * randn(grade_B, grade_B)
        B = Blade(B_vectors)

        dual_B = dual(B, C)

        # --- Verify dim, grade, and norm

        @test dim(dual_B) == test_dim
        @test grade(dual_B) == grade(C) - grade_B
        @test norm(dual_B) == norm(B)

        # --- Verify that B and dual(B) are orthogonal complements

        @test LinearAlgebra.norm(transpose(basis(dual_B)) * basis(B)) <
            10 * eps(Float64)

        # --- Verify sign(dual(B, C))

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

        # Verify dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) ≈ B
        else
            @test dual(dual_B, C) ≈ -B
        end
    end

    # ------ dim(B) == dim(C), grade(B) == grade(C)

    coefficients = randn(grade(C), grade(C))
    B_vectors = basis(C) * coefficients
    B = Blade(B_vectors)

    relative_orientation =
        sign(LinearAlgebra.det(transpose(basis(B)) * basis(C)))
    expected_result = mod(grade(B), 4) < 2 ?
        Scalar(test_dim, relative_orientation * volume(B)) :
        Scalar(test_dim, -relative_orientation * volume(B))
    @test dual(B, C) ≈ expected_result

    # ------ Error cases

    for grade_B in 2:5
        # dim(B) != dim(C)
        B = Blade(randn(test_dim + 1, grade_B))
        @test_throws DimensionMismatch dual(B, C)

        # `B` not contained in `C`
        B = Blade(randn(test_dim, grade_B))
        while LinearAlgebra.norm(rejection(basis(B), C)) < sqrt(eps(Float64))
            B = Blade(randn(test_dim, grade_B))
        end
        @test_throws ErrorException dual(B, C)
    end
end
