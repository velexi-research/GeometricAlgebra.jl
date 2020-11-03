"""
Unit tests for the wedge(x, y) function

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

@testset "wedge(B, C): B or C isa Blade" begin
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

    # --- B::Real, C::Blade
    #     B::Blade, C::Real

    # Preparations
    x = test_value_1
    B = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Blade(basis(B), volume=x * volume(B))
    @test wedge(x, B) ≈ expected_result
    @test wedge(B, x) ≈ expected_result
    @test x ∧ B ≈ expected_result
    @test B ∧ x ≈ expected_result
end

@testset "wedge(B, C): B or C isa Pseudoscalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = zero(B)
    @test wedge(B, C) == expected_result
    @test B ∧ C == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch B ∧ C

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test wedge(B, C) == expected_result
    @test wedge(C, B) == expected_result
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result

    # --- B::Pseudoscalar, C::Real
    #     B::Real, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test wedge(B, C) == expected_result
    @test wedge(C, B) == expected_result
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result

    # --- B::Pseudoscalar, v::Vector
    #     v::Vector, B::Pseudoscalar

    # dim(v) == dim(B)
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        v = Vector(rand(test_dim))

        v_wedge_B = wedge(v, B)
        expected_result = wedge(Blade(v), B)
        @test v_wedge_B ≈ expected_result
        @test v ∧ B == v_wedge_B

        B_wedge_v = wedge(B, v)
        expected_result = zero(B)
        @test B_wedge_v == expected_result
        @test B ∧ v == B_wedge_v
    end

    # dim(v) != dim(B)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    v = Vector(rand(test_dim))

    @test_throws DimensionMismatch wedge(v, B)
    @test_throws DimensionMismatch v ∧ B

    @test_throws DimensionMismatch wedge(B, v)
    @test_throws DimensionMismatch B ∧ v
end

@testset "wedge(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_wedge_C = wedge(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = B
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = B
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B

    # B::Scalar, C::Zero
    # B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = Zero()
    @test C_wedge_B === expected_result
    @test C ∧ B === C_wedge_B

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_wedge_C = wedge(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = test_value_1 * test_value_2
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B

    # B::Scalar, v::Vector
    # v::Vector, B::Scalar
    B = Scalar(test_value_1)
    v = Vector(rand(test_dim))

    # Exercise functionality and check results
    expected_result = Blade(v, volume=LinearAlgebra.norm(v) * value(B))
    @test wedge(B, v) ≈ expected_result
    @test wedge(v, B) ≈ expected_result
    @test B ∧ v ≈ expected_result
    @test v ∧ B ≈ expected_result
end

@testset "wedge(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = One()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    B_wedge_C = wedge(B, C)
    expected_result = Zero()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = Zero()
    @test C_wedge_B === expected_result
    @test C ∧ B === C_wedge_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_wedge_C = wedge(B, C)
    expected_result = C
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C == B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = C
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B == C_wedge_B
end

@testset "wedge(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    @test B + C === Zero()

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_plus_C = B + C
    expected_result = C
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = C
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result
end
