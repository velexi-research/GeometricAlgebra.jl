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

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x * B ≈ expected_result
    @test B * x ≈ expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test B * C ≈ expected_result
    @test C * B ≈ expected_result

    # --- B::Scalar, v::Vector
    #     v::Vector, B::Scalar

    # Preparations
    B = Scalar(test_value_1)
    v = Vector(rand(test_dim))

    # Exercise functionality and check results
    expected_result = Blade(v, volume=LinearAlgebra.norm(v) * value(B))
    @test B * v == expected_result
    @test v * B == expected_result
end

# --- wedge(B, C), B ∧ C

@testset "wedge(B, C): B or C isa {Scalar, Real}" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test wedge(x, B) ≈ expected_result
    @test wedge(B, x) ≈ expected_result
    @test x ∧ B ≈ expected_result
    @test B ∧ x ≈ expected_result

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test wedge(B, C) ≈ expected_result
    @test wedge(C, B) ≈ expected_result
    @test B ∧ C ≈ expected_result
    @test C ∧ B ≈ expected_result

    # --- B::Scalar, v::Vector
    #     v::Vector, B::Scalar

    # Preparations
    B = Scalar(test_value_1)
    v = Vector(rand(test_dim))

    # Exercise functionality and check results
    expected_result = Blade(v, volume=LinearAlgebra.norm(v) * value(B))
    @test wedge(B, v) ≈ expected_result
    @test wedge(v, B) ≈ expected_result
    @test B ∧ v ≈ expected_result
    @test v ∧ B ≈ expected_result
end

@testset "wedge(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- B::Pseudoscalar, C::Blade
    #     B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(test_vectors)

    expected_result = zero(B)
    @test wedge(B, C) == expected_result
    @test wedge(C, B) == expected_result
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(test_vectors)
    @test_throws DimensionMismatch B ∧ C
    @test_throws DimensionMismatch C ∧ B

    # --- B::Pseudoscalar, v::Vector
    #     v::Vector, B::Pseudoscalar

    B = Pseudoscalar(test_dim, test_value_1)
    v = Vector(rand(test_dim))

    expected_result = zero(B)
    @test wedge(B, v) == expected_result
    @test wedge(v, B) == expected_result
    @test B ∧ v == expected_result
    @test v ∧ B == expected_result

    # dim(B) != length(v)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    @test_throws DimensionMismatch B ∧ v
    @test_throws DimensionMismatch v ∧ B
end

@testset "wedge(B, C): B, C::{Blade, Vector}" begin
    # --- B::Blade, C::Blade

    B_vectors = hcat([1; 1; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vectors = hcat([0; 1; 3; 0; 0],
                     [0; 0; 0; 4; 0])
    C = Blade(C_vectors)

    # dim(B) == dim(C)
    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(hcat(B_vectors, C_vectors))
    @test B_wedge_C == B ∧ C

    # dim(B) != dim(C)
    C = Blade(rand(dim(B) + 1, 2))
    @test_throws DimensionMismatch wedge(B, C)

    # --- B::Blade, C::Vector
    #     B::Vector, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vector = [0; 1; 3; 0; 0]
    C = Blade(C_vector)

    # dim(B) == dim(C)
    expected_result = Blade(hcat(B_vectors, C_vector))

    B_wedge_C = wedge(B, C_vector)
    @test B_wedge_C ≈ expected_result
    @test B_wedge_C == B ∧ C_vector

    C_wedge_B = wedge(C_vector, B)
    @test C_wedge_B ≈ (-1)^(grade(B)) * expected_result
    @test C_wedge_B == C_vector ∧ B

    # dim(B) != dim(C)
    C = Vector(rand(dim(B) + 1))
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch wedge(C, B)

    # --- B::Vector, C::Vector

    # Preparations
    B_vector = Vector{Float64}([1, 2, 3, 0, 0])
    B = Blade(B_vector)

    C_vector = [0; 1; 3; 0; 0]
    C = Blade(C_vector)

    # dim(B) == dim(C)
    B_wedge_C = wedge(B_vector, C_vector)
    expected_result = Blade(hcat(B_vector, C_vector))
    @test B_wedge_C ≈ expected_result
    @test B_vector ∧ C_vector ≈ expected_result

    # dim(B) != dim(C)
    C_vector = Vector(rand(dim(B) + 1))
    @test_throws DimensionMismatch wedge(B_vector, C_vector)
    @test_throws DimensionMismatch wedge(C_vector, B_vector)
end

# --- proj(B, C)

@testset "proj(B, C): B or C isa Vector" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 20

    # Test vector
    v = Vector(rand(test_dim))

    # Test value
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    # ------ proj(v::Vector{<:Real}, B::Scalar)

    # return_blade == true
    B = Scalar(test_value)
    @test proj(v, B) == zero(B)

    # return_blade == false
    @test proj(v, B, return_blade=false) == 0

    # ------ proj(B::Scalar, v::Vector{<:Real})

    # return_blade == true
    B = Scalar(test_value)
    @test proj(B, v) === B

    # return_blade == false
    B = Scalar(test_value)
    @test proj(B, v, return_blade=false) == value(B)

    # ------ proj(v::Vector{<:Real}, B::Pseudoscalar)
    #        proj(B::Pseudoscalar, v::Vector{<:Real})

    # length(v) == dim(B), return_blade == true
    B = Pseudoscalar(test_dim, test_value)
    @test proj(v, B) == Blade(v)
    @test proj(B, v) == zero(B)

    # length(v) == dim(B), return_blade == false
    @test proj(v, B, return_blade=false) == v
    @test proj(B, v, return_blade=false) == 0

    # length(v) != dim(B)
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch proj(v, B)
    @test_throws DimensionMismatch proj(B, v)

    # ------ proj(v::Vector{<:Real}, B::Blade)
    #        proj(B::Blade, v::Vector{<:Real})

    # grade(B) > 1, return_blade == true
    B = Blade(rand(test_dim, 5))
    projection_vectors = basis(B) * transpose(basis(B)) * v
    @test proj(v, B) ≈ Blade(projection_vectors)
    @test proj(B, v) == zero(B)

    # grade(B) > 1, return_blade == false
    @test proj(v, B, return_blade=false) ≈ projection_vectors
    @test proj(B, v, return_blade=false) == 0

    # grade(B) == 1, return_blade == true
    B = Blade(rand(test_dim, 1))
    projection_vectors = LinearAlgebra.dot(v, basis(B)) * basis(B)
    @test proj(v, B) ≈ Blade(projection_vectors)

    projection_vectors = LinearAlgebra.dot(v, basis(B)) * v
    @test proj(B, v) ≈ Blade(projection_vectors)

    # grade(B) == 1, return_blade == false
    projection_vectors = LinearAlgebra.dot(v, basis(B)) * basis(B)
    @test proj(v, B, return_blade=false) ≈ projection_vectors

    projection_vectors = LinearAlgebra.dot(v, basis(B)) * v
    @test proj(B, v, return_blade=false) ≈ projection_vectors

    # length(v) != dim(B)
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch proj(v, B)
    @test_throws DimensionMismatch proj(B, v)
end

@testset "proj(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Scalar(value(B))
    @test proj(B, C) == expected_result

    expected_result = zero(C)
    @test proj(C, B) == expected_result
end

@testset "proj(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- B::Pseudoscalar, C::Blade
    #     B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(test_vectors)

    expected_result = zero(B)
    @test proj(B, C) == expected_result

    expected_result = C
    @test proj(C, B) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(test_vectors)
    @test_throws DimensionMismatch proj(B, C)
    @test_throws DimensionMismatch proj(C, B)
end

@testset "proj(B, C): B, C::Blade" begin
    # --- Preparations

    test_dim = 15

    # --- Valid arguments

    # ------ grade(B) < grade(C)

    B = Blade(rand(test_dim, 5))
    C = Blade(rand(test_dim, 7))

    projection = proj(B, C)
    @test projection isa Blade
    @test dim(projection) == dim(B)

    # Check norm
    norm_projection = norm(B) *
        norm(Blade(basis(C) * transpose(basis(C)) * basis(B)))
    @test norm(projection) ≈ norm_projection

    # Check that proj(B, C) is contained in C
    projection_coefficients = transpose(basis(C)) * basis(projection)
    @test LinearAlgebra.norm(projection_coefficients)^2 ≈ grade(B)

    # Check that norm(proj(B, C)
    #
    # ------ grade(B) > grade(C)

    B = Blade(rand(test_dim, 10))
    C = Blade(rand(test_dim, 3))

    @test proj(B, C) == zero(B)

    # --- Invalid arguments

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch proj(B, C)

    # --- Check consistency with proj(v::Vector, B::Blade)

    v = Vector(rand(test_dim))
    B = Blade(rand(test_dim, 3))
    @test proj(v, B) ≈ proj(Blade(v), B)
end

@testset "dual(B) tests" begin
    # --- Preparations

    # Test values
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- B::Blade

    for dim_B in 10:13
        for grade_B in 2:5
            B = Blade(rand(dim_B, grade_B))

            dual_B = dual(B)

            # --- Check dim, grade, and norm

            @test dim(dual_B) == dim_B
            @test grade(dual_B) == dim_B - grade_B
            @test norm(dual_B) == norm(B)

            # --- Check that B and dual(B) are orthogonal complements

            @test LinearAlgebra.norm(transpose(basis(dual_B)) * basis(B)) <
                10 * eps(Float64)

            # --- Check sign(dual(B))

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

            # Check dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
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

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)

    # Exercise functionality and check results
    for grade_C in 1:4
        expected_result = zero(C)
        @test dual(C, B) == expected_result
    end
end

@testset "dual(B, C): B or C isa Pseudoscalar" begin
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

    # ------ dim(B) == dim(C)

    B = Blade(rand(test_dim, 3))
    C = Pseudoscalar(test_dim, test_value_1)

    expected_result = zero(B)
    @test dual(C, B) == expected_result

    for grade_B in 2:5
        B = Blade(rand(test_dim, grade_B))
        C = Pseudoscalar(test_dim, test_value_1)

        dual_B = dual(B, C)
        @test dual_B == dual(B)

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) ≈ B
        else
            @test dual(dual_B, C) ≈ -B
        end
    end

    # ------ Error cases

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Pseudoscalar(test_dim + 1, test_value_1)
    @test_throws DimensionMismatch dual(B, C)
    @test_throws DimensionMismatch dual(C, B)
end

@testset "dual(B, C): B, C::Blade" begin
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
        while LinearAlgebra.norm(rejection(basis(B), C)) < sqrt(eps(Float64))
            B = Blade(rand(test_dim, grade_B))
        end
        @test_throws ErrorException dual(B, C)
    end
end
