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

    expected_result = Scalar(value(B) * value(C))
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result

    # --- B::Scalar, C::Blade

    # grade(C) == 1
    B = Scalar(test_value_1)
    C = Blade(test_vector)

    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(C) > 1
    B = Scalar(test_value_1)
    C = Blade(test_basis)

    expected_result = Blade(basis(C), volume=value(B) * volume(C))
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # --- B::Blade, C::Scalar

    # grade(C) == 1
    B = Blade(test_vector)
    C = Scalar(test_value_2)

    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result

    # grade(C) > 1
    B = Blade(test_basis)
    C = Scalar(test_value_2)

    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result

    # --- B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = Pseudoscalar(test_dim, value(B) * value(C))
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result

    # --- B::Pseudoscalar, C::Scalar

    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result
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

    expected_result = Scalar(x * value(B))
    @test dot(x, B) == expected_result
    @test dot(B, x) == expected_result
    @test x ⋅ B == expected_result
    @test B ⋅ x == expected_result

    # --- dot(x::Real, B::Blade)

    # grade(B) == 1
    x = test_value_1
    B = Blade(test_vector)

    expected_result = Blade(basis(B), volume=x * volume(B))
    @test dot(x, B) ≈ expected_result
    @test x ⋅ B ≈ expected_result

    # grade(B) > 1
    x = test_value_1
    B = Blade(test_basis)

    expected_result = Blade(basis(B), volume=x * volume(B))
    @test dot(x, B) ≈ expected_result
    @test x ⋅ B ≈ expected_result

    # --- dot(B::Blade, x::Real)

    # grade(B) == 1
    B = Blade(test_vector)
    x = test_value_2

    expected_result = zero(B)
    @test dot(B, x) ≈ expected_result
    @test B ⋅ x ≈ expected_result

    # grade(B) > 1
    B = Blade(test_basis)
    x = test_value_2

    expected_result = zero(B)
    @test dot(B, x) ≈ expected_result
    @test B ⋅ x ≈ expected_result

    # --- dot(x::Real, B::Pseudoscalar)

    x = test_value_1
    B = Pseudoscalar(test_dim, test_value_2)

    expected_result = Pseudoscalar(test_dim, x * value(B))
    @test dot(x, B) == expected_result
    @test x ⋅ B == expected_result

    # --- dot(B::Pseudoscalar, x::Real)

    B = Pseudoscalar(test_dim, test_value_1)
    x = test_value_2

    expected_result = zero(B)
    @test dot(B, x) == expected_result
    @test B ⋅ x == expected_result
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

    expected_result = zero(B)
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

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

        expected_result = mod(test_dim, 4) < 2 ?
            value(C) * dual(B) :
           -value(C) * dual(B)

        @test dot(B, C) ≈ expected_result
        @test B ⋅ C ≈ expected_result
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

        expected_result = mod(grade(C), 4) < 2 ?
            Scalar(value(B) * value(C)) :
            Scalar(-value(B) * value(C))

        @test dot(B, C) == expected_result
        @test B ⋅ C == expected_result
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dot(B, C)
    @test_throws DimensionMismatch B ⋅ C
end

@testset "dot(B, C): B, C::Blade" begin
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
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

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

        @test dot(B, C) ≈ expected_result
        @test B ⋅ C ≈ expected_result
    end

    # grade(B) > 1, grade(C) == 1
    B = Blade(test_basis_1)
    C = Blade(test_vector_2)

    expected_result = zero(B)
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(B) == grade(C) > 1
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))
        C = Blade(rand(test_dim, test_grade))

        expected_result = volume(B) * volume(C) *
            LinearAlgebra.det(transpose(basis(B)) * basis(C))

        @test dot(B, C) ≈ expected_result
        @test B ⋅ C ≈ expected_result
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

        @test dot(B, C) ≈ expected_result
        @test B ⋅ C ≈ expected_result
    end

    # grade(B) > grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_2)
    C = Blade(test_basis_1)

    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result

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

    # --- v::vector, B::Scalar
    #     B::Scalar, v::vector

    v = test_vector
    B = Scalar(test_value)

    expected_result = zero(B)
    @test dot(v, B) == expected_result
    @test v ⋅ B == expected_result

    expected_result = Blade(v, volume=norm(v) * test_value)
    @test dot(B, v) == expected_result
    @test B ⋅ v == expected_result

    # --- grade(B) == 1
    #     v::vector, B::Blade
    #     B::Blade, v::vector

    v = test_vector

    test_vector_B = rand(test_dim)
    B = Blade(test_vector_B)
    expected_result = LinearAlgebra.dot(test_vector, test_vector_B)
    @test dot(v, B) ≈ expected_result
    @test dot(B, v) ≈ expected_result
    @test v ⋅ B ≈ expected_result
    @test B ⋅ v ≈ expected_result

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

        expected_result = dot(Blade(v), B)
        @test dot(v, B) ≈ expected_result
        @test v ⋅ B ≈ expected_result

        expected_result = zero(B)
        @test dot(B, v) == expected_result
        @test B ⋅ v == expected_result
    end

    # length(v) != dim(B)
    v = test_vector
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch dot(v, B)
    @test_throws DimensionMismatch v ⋅ B
    @test_throws DimensionMismatch dot(B, v)
    @test_throws DimensionMismatch B ⋅ v

    # --- v::vector, B::Pseudoscalar

    # dim(v) == dim(B)
    for test_dim in 5:8
        v = Vector(rand(test_dim))
        B = Pseudoscalar(test_dim, test_value)
        expected_result = dot(Blade(v), B)
        @test v ⋅ B ≈ expected_result
        @test dot(v, B) ≈ expected_result
    end

    # dim(v) != dim(B)
    v = test_vector
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch dot(v, B)
    @test_throws DimensionMismatch v ⋅ B

    # --- B::Pseudoscalar, v::vector

    # dim(v) == dim(B)
    for test_dim in 5:8
        v = Vector(rand(test_dim))
        B = Pseudoscalar(test_dim, test_value)
        expected_result = zero(B)
        @test B ⋅ v == expected_result
        @test dot(B, v) == expected_result
    end

    # dim(v) != dim(B)
    v = test_vector
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch dot(B, v)
    @test_throws DimensionMismatch B ⋅ v
end

# --- *(B, C)
