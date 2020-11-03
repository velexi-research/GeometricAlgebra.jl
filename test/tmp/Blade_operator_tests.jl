"""
Operator unit tests for Blade type.

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

# ------ Unary operations

@testset "Blade: -(B)" begin
    # Preparations
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)

    # B::Blade
    negative_B = Blade(B, volume=-volume(B))
    @test -B == negative_B

    @test -negative_B == B
end

@testset "Blade: reverse(B)" begin
    # mod(grade, 4) == 1
    vectors = Vector([3; 4; 0; 0; 0])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 2
    vectors = Matrix([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2
end

@testset "Blade: dual(B)" begin
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

@testset "Blade: reciprocal(B)" begin
    # mod(grade, 4) == 1
    vectors = Vector([3; 4; 0; 0; 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 2
    vectors = Matrix([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1
end

# ------ Binary operations

@testset "Blade: proj(B, C)" begin
    # --- Preparations

    test_dim = 15

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vectors = rand(test_dim, 3)

    # --- B::Blade, C::Blade

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

    # Check that norm(proj(B, C)) TODO

    # ------ grade(B) > grade(C)

    B = Blade(rand(test_dim, 10))
    C = Blade(rand(test_dim, 3))

    @test proj(B, C) == zero(B)

    # ------ dim(B) != dim(C)

    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch proj(B, C)

    # ------ Check consistency with proj(v::Vector, B::Blade)

    v = Vector(rand(test_dim))
    B = Blade(rand(test_dim, 3))
    @test proj(v, B) ≈ proj(Blade(v), B)

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

@testset "Blade: dual(B, C)" begin
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

@testset "Blade: dual(B, C)" begin
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

    C = Scalar(test_value_1)

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

    C = One()

    for test_grade in 5:8
        B = Blade(randn(test_dim, test_grade))

        B_dual_C = dual(B, C)
        expected_result = mod(test_grade, 4) < 2 ?
            Blade(C, volume=1, copy_basis=false) :
            Blade(C, volume=-1, copy_basis=false)
        @test B_dual_C isa Blade
        @test B_dual_C == expected_result

        C_dual_B = dual(C, B)
        expected_result = Zero(C)
        @test C_dual_B === expected_result
    end
end

#------------------------- NEEDS REORGANIZATION ---------------------#

# --- dot(B, C), B ⋅ C

@testset "dot(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vector = rand(test_dim)

    # Test basis
    test_grade = 3
    test_basis = rand(test_dim, test_grade)

    # --- B::Scalar, C::Scalar

    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_dot_C = dot(B, C)
    expected_result = Scalar(value(B) * value(C))
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # --- B::Scalar, C::Blade

    # grade(C) == 1
    B = Scalar(test_value_1)
    C = Blade(test_vector)

    B_dot_C = dot(B, C)
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test B_dot_C ≈ expected_result
    @test B ⋅ C == B_dot_C

    # grade(C) > 1
    B = Scalar(test_value_1)
    C = Blade(test_basis)

    B_dot_C = dot(B, C)
    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C == B_dot_C

    # --- B::Blade, C::Scalar

    # grade(C) == 1
    B = Blade(test_vector)
    C = Scalar(test_value_2)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # grade(C) > 1
    B = Blade(test_basis)
    C = Scalar(test_value_2)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # --- B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_dot_C = dot(B, C)
    expected_result = Pseudoscalar(test_dim, value(B) * value(C))
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # --- B::Pseudoscalar, C::Scalar

    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C
end

@testset "dot(B, C): B or C isa Real" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vector = rand(test_dim)

    # Test basis
    test_grade = 3
    test_basis = rand(test_dim, test_grade)

    # --- dot(x::Real, B::Scalar)
    #     dot(B::Scalar, x::Real)

    x = test_value_1
    B = Scalar(test_value_2)

    x_dot_B = dot(x, B)
    expected_result = Scalar(x * value(B))
    @test x_dot_B == expected_result
    @test x ⋅ B == x_dot_B

    B_dot_x = dot(B, x)
    expected_result = Scalar(x * value(B))
    @test B_dot_x == expected_result
    @test B ⋅ x == B_dot_x

    # --- dot(x::Real, B::Blade)

    # grade(B) == 1
    x = test_value_1
    B = Blade(test_vector)

    x_dot_B = dot(x, B)
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x_dot_B ≈ expected_result
    @test x ⋅ B == x_dot_B

    # grade(B) > 1
    x = test_value_1
    B = Blade(test_basis)

    x_dot_B = dot(x, B)
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test x_dot_B ≈ expected_result
    @test x ⋅ B == x_dot_B

    # --- dot(B::Blade, x::Real)

    # grade(B) == 1
    B = Blade(test_vector)
    x = test_value_2

    B_dot_x = dot(B, x)
    expected_result = zero(B)
    @test B_dot_x == expected_result
    @test B ⋅ x == B_dot_x

    # grade(B) > 1
    B = Blade(test_basis)
    x = test_value_2

    B_dot_x = dot(B, x)
    expected_result = zero(B)
    @test B_dot_x == expected_result
    @test B ⋅ x == B_dot_x

    # --- dot(x::Real, B::Pseudoscalar)

    x = test_value_1
    B = Pseudoscalar(test_dim, test_value_2)

    x_dot_B = dot(x, B)
    expected_result = Pseudoscalar(test_dim, x * value(B))
    @test x_dot_B == expected_result
    @test x ⋅ B == x_dot_B

    # --- dot(B::Pseudoscalar, x::Real)

    B = Pseudoscalar(test_dim, test_value_1)
    x = test_value_2

    B_dot_x = dot(B, x)
    expected_result = zero(B)
    @test B_dot_x == expected_result
    @test B ⋅ x == B_dot_x
end

@testset "dot(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Pseudoscalar, C::Blade

    # dim(B) == dim(C)
    test_dim = 10
    test_grade = 5
    test_basis = rand(test_dim, test_grade)

    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(test_basis)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(test_basis)
    @test_throws DimensionMismatch dot(B, C)
    @test_throws DimensionMismatch B ⋅ C

    # --- B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    test_grade = 3
    for test_dim in 10:13
        B = Blade(rand(test_dim, test_grade))
        C = Pseudoscalar(test_dim, test_value_2)

        B_dot_C = dot(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            value(C) * dual(B) :
           -value(C) * dual(B)

        @test B_dot_C ≈ expected_result
        @test B ⋅ C == B_dot_C
    end

    # dim(B) != dim(C)
    B = Blade(test_basis)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dot(B, C)
    @test_throws DimensionMismatch B ⋅ C

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    for test_dim in 10:13
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_dot_C = dot(B, C)

        expected_result = mod(grade(C), 4) < 2 ?
            Scalar(value(B) * value(C)) :
            Scalar(-value(B) * value(C))

        @test B_dot_C == expected_result
        @test B ⋅ C == B_dot_C
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dot(B, C)
    @test_throws DimensionMismatch B ⋅ C
end

@testset "Blade: contractl(B, C)" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10
    test_grade_1 = 3
    test_grade_2 = 5

    # Test vectors
    test_vector_1 = rand(test_dim)
    test_vector_2 = rand(test_dim)

    # Test bases
    test_basis_1 = rand(test_dim, test_grade_1)
    test_basis_2 = rand(test_dim, test_grade_2)

    # --- Test cases

    # grade(B) == grade(C) == 1
    B = Blade(test_vector_1)
    C = Blade(test_vector_2)

    expected_result =
        Scalar(LinearAlgebra.dot(basis(B), basis(C)) * volume(B) * volume(C))

    B_dot_C = dot(B, C)
    @test B_dot_C ≈ expected_result
    @test B ⋅ C == B_dot_C

    # grade(B) == 1, grade(C) > 1
    B = Blade(test_vector_1)

    for test_grade in 5:8
        C = Blade(test_basis_2)

        F = LinearAlgebra.qr(test_basis_2)
        Q = Matrix(F.Q)
        projection = Blade(Q * transpose(Q) * test_vector_1)
        expected_volume_C = prod(LinearAlgebra.diag(F.R))
        expected_result = mod(grade(C), 4) < 2 ?
            expected_volume_C * dual(projection, C) :
           -expected_volume_C * dual(projection, C)

        B_dot_C = dot(B, C)
        @test B_dot_C ≈ expected_result
        @test B ⋅ C == B_dot_C
    end

    # grade(B) > 1, grade(C) == 1
    B = Blade(test_basis_1)
    C = Blade(test_vector_2)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == B_dot_C

    # grade(B) == grade(C) > 1
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))
        C = Blade(rand(test_dim, test_grade))

        expected_result = volume(B) * volume(C) *
            LinearAlgebra.det(transpose(basis(B)) * basis(C))

        B_dot_C = dot(B, C)
        @test B_dot_C ≈ expected_result
        @test B ⋅ C == B_dot_C
    end

    # grade(B) < grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_1)

    for test_grade in 5:8
        test_basis_C = rand(test_dim, test_grade)
        C = Blade(test_basis_C)

        F = LinearAlgebra.qr(test_basis_C)
        Q = Matrix(F.Q)
        projection = Blade(Q * transpose(Q) * test_basis_1)
        expected_volume_C = prod(LinearAlgebra.diag(F.R))
        expected_result = mod(grade(C), 4) < 2 ?
            expected_volume_C * dual(projection, C) :
           -expected_volume_C * dual(projection, C)

        B_dot_C = dot(B, C)
        @test B_dot_C ≈ expected_result
        @test B ⋅ C == B_dot_C
    end

    # grade(B) > grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_2)
    C = Blade(test_basis_1)

    B_dot_C = dot(B, C)
    expected_result = zero(B)
    @test B_dot_C == expected_result
    @test B ⋅ C == B_dot_C

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch dot(B, C)
    @test_throws DimensionMismatch B ⋅ C
end

@testset "dot(B, C): B or C isa Vector" begin
    # Notes: these tests are based on the assumption that dot(B, C) is
    #        correct when B and C are both AbstractBlades.

    # --- Preparations

    # Dimension of embedding space
    test_dim = 10
    test_grade = 3

    # Test value
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # Test vectors
    test_vector = rand(test_dim)

    # --- grade(B) == 1
    #     v::vector, B::Blade
    #     B::Blade, v::vector

    v = test_vector

    test_vector_B = rand(test_dim)
    B = Blade(test_vector_B)

    expected_result = LinearAlgebra.dot(test_vector, test_vector_B)

    B_dot_v = dot(B, v)
    @test B_dot_v ≈ expected_result
    @test B ⋅ v == B_dot_v

    v_dot_B = dot(v, B)
    @test v_dot_B ≈ expected_result
    @test v ⋅ B == v_dot_B

    # length(v) != dim(B)
    v = test_vector
    B = Blade(rand(test_dim + 1))
    @test_throws DimensionMismatch dot(v, B)
    @test_throws DimensionMismatch v ⋅ B
    @test_throws DimensionMismatch dot(B, v)
    @test_throws DimensionMismatch B ⋅ v

    # --- grade(B) > 1
    #     v::vector, B::Blade
    #     B::Blade, v::vector

    v = test_vector
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))

        v_dot_B = dot(v, B)
        expected_result = dot(Blade(v), B)
        @test v_dot_B ≈ expected_result
        @test v ⋅ B == v_dot_B

        B_dot_v = dot(B, v)
        expected_result = zero(B)
        @test B_dot_v == expected_result
        @test B ⋅ v == B_dot_v
    end

    # length(v) != dim(B)
    v = test_vector
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch dot(v, B)
    @test_throws DimensionMismatch v ⋅ B
    @test_throws DimensionMismatch dot(B, v)
    @test_throws DimensionMismatch B ⋅ v
end

# ------ Comparison operations

@testset "Blade: ==(B, C)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) == volume(C)
    # basis(B) == basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            if precision_type1 == precision_type2
                @test B == C
            else
                @test B != C
            end
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # volume(B) != volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=2*volume(B))
            @test B != C
        end
    end

    # basis(B) != basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 5; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                      volume=test_volume)
            @test B != C
        end
    end
end

@testset "Blade: ≈(B, C)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) ≈ volume(C)
    # basis(B) ≈ basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            @test B ≈ C
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # volume(B) ≉ volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=2*volume(B))
            @test B ≉ C
        end
    end

    # basis(B) ≉ basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 10; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                      volume=test_volume)
            @test B ≉ C
        end
    end

    # B and C have opposite orientations
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=-test_volume)
            @test B ≉ C
        end
    end
end
