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
    C = Blade(test_basis_2)

    F = LinearAlgebra.qr(test_basis_2)
    Q = Matrix(F.Q)
    projection = Blade(Q * transpose(Q) * test_vector_1)
    expected_volume_C = prod(LinearAlgebra.diag(F.R))
    expected_result = expected_volume_C * dual(projection, C)
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(B) > 1, grade(C) == 1
    B = Blade(test_basis_1)
    C = Blade(test_vector_2)

    expected_result = zero(B)
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(B) == grade(C) > 1
    B = Blade(test_basis_1)
    C = Blade(rand(test_dim, test_grade_1))

    expected_result =
        volume(B) * volume(C) * volume(Blade(transpose(basis(B)) * basis(C)))
    expected_result = mod(grade(C), 4) < 2 ?
        expected_result : -expected_result

    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(B) < grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_1)
    C = Blade(test_basis_2)

    F = LinearAlgebra.qr(test_basis_2)
    Q = Matrix(F.Q)
    projection = Blade(Q * transpose(Q) * test_basis_1)
    expected_volume_C = prod(LinearAlgebra.diag(F.R))
    expected_result = expected_volume_C * dual(projection, C)
    @test dot(B, C) ≈ expected_result
    @test B ⋅ C ≈ expected_result

    # grade(B) > grade(C), grade(B) > 1, grade(C) > 1,
    B = Blade(test_basis_2)
    C = Blade(test_basis_1)

    expected_result = zero(B)
    @test dot(B, C) == expected_result
    @test B ⋅ C == expected_result
end

@testset "dot(B, C): B or C isa Vector" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test vectors
    v = rand(test_dim)
    w = rand(test_dim)

    # --- v::vector, B::Blade
    #     B::Blade, v::vector

    # Preparations
    B = Blade(w)

    # Exercise functionality and check results
    expected_result = LinearAlgebra.dot(v, w)
    @test v ⋅ B ≈ expected_result
    @test B ⋅ v ≈ expected_result
    @test dot(v, B) ≈ expected_result
    @test dot(B, v) ≈ expected_result

    # --- v::vector, B::Scalar
    #     B::Scalar, v::vector

    # TODO
end

@testset "dot(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Dimension of embedding space
    test_dim = 10

    # Test vectors
    v = rand(test_dim)
    w = rand(test_dim)

    # --- v::vector, B::Pseudoscalar
    #     B::Pseudoscalar, v::vector

    # TODO
end

# --- *(B, C)
