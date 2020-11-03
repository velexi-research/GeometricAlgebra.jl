"""
Unit tests for the contractl(B, C) function

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

@testset "contractl(B, C): B or C isa Blade" begin
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

    # --- B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    test_grade = 3
    for test_dim in 10:13
        B = Blade(rand(test_dim, test_grade))
        C = Pseudoscalar(test_dim, test_value_2)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            value(C) * dual(B) :
           -value(C) * dual(B)

        @test B_contractl_C ≈ expected_result
        @test B < C == B_contractl_C
    end

    # dim(B) != dim(C)
    B = Blade(test_basis)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C

    # --- B::Pseudoscalar, C::Blade

    # dim(B) == dim(C)
    test_dim = 10
    test_grade = 5
    test_basis = rand(test_dim, test_grade)

    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(test_basis)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(test_basis)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C

    # --- B::Blade, C::Scalar
    #     B::Scalar, B::Blade

    # grade(C) == 1
    B = Blade(test_vector)
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Blade(basis(B), volume=value(C) * volume(B))
    @test C_contractl_B ≈ expected_result
    @test C < B == C_contractl_B

    # grade(C) > 1
    B = Blade(test_basis)
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Blade(basis(B), volume=value(C) * volume(B))
    @test C_contractl_B ≈ expected_result
    @test C < B == C_contractl_B

    # --- B::Blade, C::Real
    #     C::Real, B::Blade

    # grade(B) == 1
    B = Blade(test_vector)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Blade(basis(B), volume=C * volume(B))
    @test C_contractl_B ≈ expected_result
    @test C < B == C_contractl_B

    # grade(B) > 1
    B = Blade(test_basis)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Blade(basis(B), volume=C * volume(B))
    @test C_contractl_B ≈ expected_result
    @test C < B == C_contractl_B

    # --- v::Vector, B::Blade
    #     B::Blade, v::Vector

    # grade(B) == 1
    v = test_vector

    test_vector_B = rand(test_dim)
    B = Blade(test_vector_B)

    expected_result = LinearAlgebra.contractl(test_vector, test_vector_B)

    B_contractl_v = contractl(B, v)
    @test B_contractl_v ≈ expected_result
    @test B < v == B_contractl_v

    v_contractl_B = contractl(v, B)
    @test v_contractl_B ≈ expected_result
    @test v < B == v_contractl_B

    # length(v) != dim(B)
    v = test_vector
    B = Blade(rand(test_dim + 1))
    @test_throws DimensionMismatch contractl(v, B)
    @test_throws DimensionMismatch v < B
    @test_throws DimensionMismatch contractl(B, v)
    @test_throws DimensionMismatch B < v

    # grade(B) > 1
    v = test_vector
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))

        v_contractl_B = contractl(v, B)
        expected_result = contractl(Blade(v), B)
        @test v_contractl_B ≈ expected_result
        @test v < B == v_contractl_B

        B_contractl_v = contractl(B, v)
        expected_result = zero(B)
        @test B_contractl_v == expected_result
        @test B < v == B_contractl_v
    end

    # length(v) != dim(B)
    v = test_vector
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch contractl(v, B)
    @test_throws DimensionMismatch v < B
    @test_throws DimensionMismatch contractl(B, v)
    @test_throws DimensionMismatch B < v
end

@testset "contractl(B, C): B::Blade, C::Blade" begin
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
        Scalar(LinearAlgebra.contractl(basis(B), basis(C)) *
               volume(B) * volume(C))

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ expected_result
    @test B < C == B_contractl_C

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

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ expected_result
        @test B < C == B_contractl_C
    end

    # grade(B) > 1, grade(C) == 1
    B = Blade(test_basis_1)
    C = Blade(test_vector_2)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test contractl(B, C) == expected_result
    @test B < C == B_contractl_C

    # grade(B) == grade(C) > 1
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))
        C = Blade(rand(test_dim, test_grade))

        expected_result = volume(B) * volume(C) *
            LinearAlgebra.det(transpose(basis(B)) * basis(C))

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ expected_result
        @test B < C == B_contractl_C
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

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ expected_result
        @test B < C == B_contractl_C
    end

    # grade(B) > grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_2)
    C = Blade(test_basis_1)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

@testset "contractl(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    for test_dim = 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_contractl_C isa Scalar
        @test B_contractl_C == expected_result
        @test (B < C) == B_contractl_C
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C

    # --- B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = Pseudoscalar(test_dim, value(B) * value(C))
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    # --- B::Pseudoscalar, C::Scalar

    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    # --- B::Pseudoscalar, C::Real
    #     C::Pseudoscalar, B::Real

    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = zero(B)
    @test B_contractl_C == expected_result
    @test B < C == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Pseudoscalar(test_dim, C * value(B))
    @test C_contractl_B == expected_result
    @test C < B == C_contractl_B

    # --- B::Pseudoscalar, v::Vector
    #     v::Vector, B::Pseudoscalar

    # dim(v) == dim(B)
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        v = Vector(rand(test_dim))

        v_contractl_B = contractl(v, B)
        expected_result = contractl(Blade(v), B)
        @test v_contractl_B ≈ expected_result
        @test v < B == v_contractl_B

        B_contractl_v = contractl(B, v)
        expected_result = zero(B)
        @test B_contractl_v == expected_result
        @test B < v == B_contractl_v
    end

    # dim(v) != dim(B)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    v = test_vector

    @test_throws DimensionMismatch contractl(v, B)
    @test_throws DimensionMismatch v < B

    @test_throws DimensionMismatch contractl(B, v)
    @test_throws DimensionMismatch B < v
end

@testset "contractl(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    # --- B::Scalar, C::One
    #     B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_contractl_C = contractl(B, C)
    expected_result = B
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = B
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B

    # --- B::Scalar, C::Zero
    #     B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B

    # --- B::Scalar, C::Real
    #     B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B

    # --- B::Scalar, v::Vector
    #     v::Vector, B::Scalar

    B = Scalar(test_value_1)
    v = test_vector

    v_contractl_B = contractl(v, B)
    expected_result = zero(B)
    @test v_contractl_B == expected_result
    @test (v < B) == v_contractl_B

    B_contractl_v = contractl(B, v)
    expected_result = Blade(v, volume=norm(v) * test_value)
    @test B_contractl_v == expected_result
    @test (B < v) == B_contractl_v
end

@testset "contractl(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_contractl_C = contractl(B, C)
    expected_result = One()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_contractl_C = contractl(B, C)
    expected_result = C
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) == B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = C
    @test C_contractl_B isa Scalar
    @test C_contractl_B == expected_result
    @test (C < B) == C_contractl_B
end

@testset "contractl(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_contractl_C = contractl(B, C)
    expected_result = Zero()
    @test B_contractl_C === expected_result
    @test (B < C) === B_contractl_C

    C_contractl_B = contractl(C, B)
    expected_result = Zero()
    @test C_contractl_B === expected_result
    @test (C < B) === C_contractl_B
end
