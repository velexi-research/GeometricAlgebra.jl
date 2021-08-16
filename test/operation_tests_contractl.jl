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
using LinearAlgebra: ⋅, det, diag, qr
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- File inclusions

# Test utilities
include("test_utils.jl")

#=
# --- Tests

# ------ M::Multivector

@testset "contractl(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "contractl(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "contractl(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Blade, C::Blade)" begin
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

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Scalar((basis(B) ⋅ basis(C)) * volume(B) * volume(C))
    @test (B < C) == B_contractl_C

    # grade(B) == 1, grade(C) > 1
    B = Blade(test_vector_1)

    for test_grade in 5:8
        C = Blade(test_basis_2)

        F = qr(test_basis_2)
        Q = Matrix(F.Q)
        projection = Blade(Q * transpose(Q) * test_vector_1)
        expected_volume_C = prod(diag(F.R))
        expected_result = mod(grade(C), 4) < 2 ?
            expected_volume_C * dual(projection, C) :
           -expected_volume_C * dual(projection, C)

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ expected_result
        @test (B < C) == B_contractl_C
    end

    # grade(B) > 1, grade(C) == 1
    B = Blade(test_basis_1)
    C = Blade(test_vector_2)

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) == grade(C) > 1
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))
        C = Blade(rand(test_dim, test_grade))

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ volume(B) * volume(C) *
                              det(transpose(basis(B)) * basis(C))

        @test (B < C) == B_contractl_C
    end

    # grade(B) < grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_1)

    for test_grade in 5:8
        test_basis_C = rand(test_dim, test_grade)
        C = Blade(test_basis_C)

        F = qr(test_basis_C)
        Q = Matrix(F.Q)
        projection = Blade(Q * transpose(Q) * test_basis_1)
        expected_volume_C = prod(diag(F.R))
        expected_result = mod(grade(C), 4) < 2 ?
            expected_volume_C * dual(projection, C) :
           -expected_volume_C * dual(projection, C)

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ expected_result
        @test (B < C) == B_contractl_C
    end

    # grade(B) > grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_2)
    C = Blade(test_basis_1)

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # dim(B) != dim(C)
    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

@testset "contractl(B::Blade, C::Pseudoscalar)" begin
    # ---Preparations

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # dim(B) == dim(C)
    test_grade = 3
    for test_dim in 10:13
        B = Blade(rand(test_dim, test_grade))
        C = Pseudoscalar(test_dim, test_value)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            test_value * dual(B) :
           -test_value * dual(B)

        @test B_contractl_C ≈ expected_result
        @test (B < C) == B_contractl_C
    end

    # dim(B) != dim(C)
    test_dim = 8
    B = Blade(rand(test_dim, test_grade))
    C = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

@testset "contractl(B::Blade, C::Scalar)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    test_vector = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1
    B = Blade(test_vector)
    C = Scalar(test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) > 1
    B = Blade(test_basis)
    C = Scalar(test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Blade, C::One)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3
    test_vector = rand(test_dim)
    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1
    B = Blade(test_vector)
    C = One()

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) > 1
    B = Blade(test_basis)
    C = One()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Blade, C::Zero)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3
    test_vector = rand(test_dim)
    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1
    B = Blade(test_vector)
    C = Zero()

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) > 1
    B = Blade(test_basis)
    C = Zero()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Blade, C::Real)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    test_vector = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests
    #     C::Real, B::Blade

    # grade(B) == 1
    B = Blade(test_vector)
    C = test_value

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) > 1
    B = Blade(test_basis)
    C = test_value

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Blade, C::Vector)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    test_vector_1 = rand(test_dim)
    test_vector_2 = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1, dim(B) == length(C)
    B = Blade(test_vector_1)
    C = test_vector_2

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ test_vector_1 ⋅ test_vector_2
    @test (B < C) == B_contractl_C

    # grade(B) == 1, dim(B) != length(C)
    B = Blade(rand(test_dim + 1))
    C = test_vector_2
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C

    # grade(B) > 1, dim(B) == length(C)
    C = test_vector_2
    for test_grade in 5:8
        B = Blade(rand(test_dim, test_grade))

        @test iszero(contractl(B, C))
        @test iszero(B < C)
    end

    # grade(B) > 1, dim(B) != length(C)
    B = Blade(rand(test_dim + 1, 3))
    C = test_vector_2
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

# ------ B::Pseudoscalar

@testset "contractl(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Pseudoscalar, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 5

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    C = Blade(rand(test_dim, test_grade))

    # --- Tests

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch (B < C)
end
=#

@testset "contractl(B::Pseudoscalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = get_random_value(1)  # add 2 to keep value away from 0
    test_value_2 = get_random_value(1)  # add 2 to keep value away from 0

    # --- Tests

    # dim(B) == dim(C)
    for test_dim = 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_contractl_C isa AbstractScalar
        @test B_contractl_C == expected_result
        @test (B < C) == B_contractl_C
    end

    # dim(B) != dim(C)
    test_dim = 5
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

#=
@testset "contractl(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 11
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Pseudoscalar, C::One)" begin
    test_dim = 12
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = One()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Pseudoscalar, C::Zero)" begin
    test_dim = 12
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = Zero()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Pseudoscalar, C::Real)" begin
    test_dim = 12
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Pseudoscalar, C::Vector)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # dim(v) == dim(B)
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value)
        C = Vector(rand(test_dim))

        @test iszero(contractl(B, C))
        @test iszero(B < C)
    end

    # dim(B) != dim(C)
    test_dim = 10
    B = Pseudoscalar(test_dim + 1, test_value)
    C = rand(test_dim)

    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

# ------ B::Scalar

@testset "contractl(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Scalar, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    test_vector = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(C) == 1
    B = Scalar(test_value)
    C = Blade(test_vector)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(basis(C), volume=test_value * volume(C))
    @test (B < C) == B_contractl_C

    # grade(C) > 1
    B = Scalar(test_value)
    C = Blade(test_basis)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(basis(C), volume=test_value * volume(C))
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Scalar, C::Pseudoscalar)" begin
    test_dim = 10

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test (B < C) == B_contractl_C
end
=#

@testset "contractl(B::Scalar, C::Scalar)" begin
    test_value_1 = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value_1 * test_value_2
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Scalar, C::One)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = One()

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == B
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Scalar, C::Zero)" begin
    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    B = Scalar(test_value)

    C = Zero()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

#=
@testset "contractl(B::Scalar, C::Real)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value_1 * test_value_2
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Scalar, C::Vector)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = rand(5)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(C, volume=norm(C) * test_value)
    @test (B < C) == B_contractl_C
end

# ------ B::One

@testset "contractl(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::One, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3
    test_vector = rand(test_dim)
    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1
    B = One()
    C = Blade(test_vector)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C == C
    @test (B < C) == B_contractl_C

    # grade(B) > 1
    B = One()
    C = Blade(test_basis)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C == C
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::One, C::Pseudoscalar)" begin
    B = One()

    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Pseudoscalar(test_dim, test_value)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C == Pseudoscalar(test_dim, test_value)
    @test (B < C) == B_contractl_C
end
=#

@testset "contractl(B::One, C::Scalar)" begin
    B = One()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::One, C::One)" begin
    B = One()
    C = One()
    @test isone(contractl(B, C))
    @test isone(B < C)
end

@testset "contractl(B::One, C::Zero)" begin
    B = One()
    C = Zero()
    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

#=
@testset "contractl(B::One, C::Real)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::One, C::Vector)" begin
    B = One()
    C = rand(5)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(C, volume=norm(C))
    @test (B < C) == B_contractl_C
end

# ------ B::Zero

@testset "contractl(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Zero, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3
    test_vector = rand(test_dim)
    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(B) == 1
    B = Zero()
    C = Blade(test_vector)

    @test iszero(contractl(B, C))
    @test iszero(B < C)

    # grade(B) > 1
    B = Zero()
    C = Blade(test_basis)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Zero, C::Pseudoscalar)" begin
    B = Zero()

    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Pseudoscalar(test_dim, test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end
=#

@testset "contractl(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    C = Scalar(test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()
    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

#=
@testset "contractl(B::Zero, C::Real)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Zero, C::Vector)" begin
    B = Zero()
    C = rand(5)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

# ------ B::Real

@testset "contractl(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Real, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    # Test values
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    test_vector = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(C) == 1
    B = test_value
    C = Blade(test_vector)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(basis(C), volume=test_value * volume(C))
    @test (B < C) == B_contractl_C

    # grade(C) > 1
    B = test_value
    C = Blade(test_basis)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ Blade(basis(C), volume=test_value * volume(C))
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Real, C::Pseudoscalar)" begin
    test_dim = 12
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value_1 * test_value_2
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_contractl_C = contractl(B, C)
    @test B_contractl_C isa AbstractScalar
    @test B_contractl_C == test_value
    @test (B < C) == B_contractl_C
end

@testset "contractl(B::Real, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

# ------ B::Vector

@testset "contractl(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "contractl(B::Vector, C::Blade)" begin
    # --- Preparations

    test_dim = 10
    test_grade = 3

    test_vector_1 = rand(test_dim)
    test_vector_2 = rand(test_dim)

    test_basis = rand(test_dim, test_grade)

    # --- Tests

    # grade(C) == 1, dim(B) == length(C)
    B = test_vector_1
    C = Blade(test_vector_2)

    B_contractl_C = contractl(B, C)
    @test B_contractl_C ≈ test_vector_1 ⋅ test_vector_2
    @test (B < C) == B_contractl_C

    # grade(C) == 1, dim(B) != length(C)
    C = Blade(rand(test_dim + 1))
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C

    # grade(C) > 1, dim(B) == length(C)
    B = test_vector_1
    for test_grade in 5:8
        C = Blade(rand(test_dim, test_grade))

        B_contractl_C = contractl(B, C)
        @test B_contractl_C ≈ contractl(Blade(B), C)
        @test (B < C) == B_contractl_C
    end

    # grade(C) > 1, dim(B) != length(C)
    B = test_vector_1
    C = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

@testset "contractl(B::Vector, C::Pseudoscalar)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # dim(B) == dim(C)
    for test_dim in 5:8
        B = Vector(rand(test_dim))
        C = Pseudoscalar(test_dim, test_value)

        B_contract_C = contractl(B, C)
        @test B_contract_C ≈ contractl(Blade(B), C)
        @test (B< C) == B_contract_C
    end

    # dim(B) != dim(C)
    test_dim = 10
    B = Pseudoscalar(test_dim + 1, test_value)
    C = rand(test_dim)

    @test_throws DimensionMismatch contractl(B, C)
    @test_throws DimensionMismatch B < C
end

@testset "contractl(B::Vector, C::Scalar)" begin
    B = rand(5)

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Vector, C::One)" begin
    B = rand(5)
    C = One()
    @test iszero(contractl(B, C))
    @test iszero(B < C)
end

@testset "contractl(B::Vector, C::Zero)" begin
    B = rand(5)
    C = Zero()
    @test iszero(contractl(B, C))
    @test iszero(B < C)
end
=#
