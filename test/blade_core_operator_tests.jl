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
using Test

# GeometricAlgebra.jl
using GeometricAlgebra


# --- *(B, C): scalar multiplication

@testset "*(B, C) tests: B or C isa {Scalar, Real}" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Blade vectors
    vectors = randn(6, 3)

    # Pseudoscalar dimension
    dim = 100

    # --- x::Real, B::Scalar
    #     B::Scalar, x::Real

    # Preparations
    x = test_value_1
    B = Scalar(test_value_2)

    # Exercise functionality and check results
    expected_x_times_B = Scalar(x * value(B))
    @test x * B == expected_x_times_B
    @test B * x == expected_x_times_B

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    # Exercise functionality and check results
    expected_B_times_C = Scalar(value(B) * value(C))
    @test B * C == expected_B_times_C

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(vectors)

    # Exercise functionality and check results
    expected_x_times_B = Blade(basis(B), volume=x * volume(B))
    @test x * B ≈ expected_x_times_B
    @test B * x ≈ expected_x_times_B

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(vectors)

    # Exercise functionality and check results
    expected_B_times_C = Blade(basis(C), volume=value(B) * volume(C))
    @test B * C ≈ expected_B_times_C
    @test C * B ≈ expected_B_times_C

    # --- x::Real, B::Pseudoscalar
    #     B::Pseudoscalar, x::Real

    # Preparations
    x = test_value_1
    B = Pseudoscalar(dim, test_value_2)

    # Exercise functionality and check results
    expected_x_times_B = Pseudoscalar(dim, x * value(B))
    @test x * B == expected_x_times_B
    @test B * x == expected_x_times_B

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Pseudoscalar(dim, test_value_2)

    # Exercise functionality and check results
    expected_B_times_C = Pseudoscalar(dim, value(B) * value(C))
    @test B * C ≈ expected_B_times_C
    @test C * B ≈ expected_B_times_C
end

# --- ∧(B, C)

@testset "∧(B, C) tests: B or C isa {Scalar, Real}" begin
    # --- Preparations

    # Test values
    test_value_1 = rand()
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand()
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Blade vectors
    vectors = randn(6, 3)

    # Pseudoscalar dimension
    dim = 100

    # --- x::Real, B::Scalar
    #     B::Scalar, x::Real

    # Preparations
    x = test_value_1
    B = Scalar(test_value_2)

    # Exercise functionality and check results
    expected_x_wedge_B = Scalar(x * value(B))
    @test x ∧ B == expected_x_wedge_B
    @test B ∧ x == expected_x_wedge_B

    # --- B::Scalar, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    # Exercise functionality and check results
    expected_B_wedge_C = Scalar(value(B) * value(C))
    @test B ∧ C == expected_B_wedge_C

    # --- x::Real, B::Blade
    #     B::Blade, x::Real

    # Preparations
    x = test_value_1
    B = Blade(vectors)

    # Exercise functionality and check results
    expected_x_wedge_B = Blade(basis(B), volume=x * volume(B))
    @test x ∧ B ≈ expected_x_wedge_B
    @test B ∧ x ≈ expected_x_wedge_B

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(vectors)

    # Exercise functionality and check results
    expected_B_wedge_C = Blade(basis(C), volume=value(B) * volume(C))
    @test B ∧ C ≈ expected_B_wedge_C
    @test C ∧ B ≈ expected_B_wedge_C

    # --- x::Real, B::Pseudoscalar
    #     B::Pseudoscalar, x::Real

    # Preparations
    x = test_value_1
    B = Pseudoscalar(dim, test_value_2)

    # Exercise functionality and check results
    expected_x_wedge_B = Pseudoscalar(dim, x * value(B))
    @test x ∧ B == expected_x_wedge_B
    @test B ∧ x == expected_x_wedge_B

    # --- B::Scalar, C::Pseudoscalar
    #     B::Pseudoscalar, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Pseudoscalar(dim, test_value_2)

    # Exercise functionality and check results
    expected_B_wedge_C = Pseudoscalar(dim, value(B) * value(C))
    @test B ∧ C ≈ expected_B_wedge_C
    @test C ∧ B ≈ expected_B_wedge_C
end

@testset "∧(B, C) tests: B, C::{Blade, Vector}" begin
    # --- B, C::Blade

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

    expected_B_wedge_C = Blade(hcat(B_vectors, C_vector))

    # Exercise functionality
    B_wedge_C = B ∧ C_vector
    @test B_wedge_C ≈ expected_B_wedge_C
    @test B_wedge_C == outer(B, C)

    C_wedge_B = C_vector ∧ B
    @test C_wedge_B ≈ (-1)^(grade(B)) * expected_B_wedge_C
    @test C_wedge_B == outer(C, B)

    # --- B::Vector, C::Vector

    # Preparations
    B_vector = Vector{Float64}([0, 2, 0, 0, 0])
    B = Blade(B_vector)

    C_vector = [0; 0; 3; 0; 0]
    C = Blade(C_vector)

    expected_B_wedge_C = Blade(hcat(B_vector, C_vector))

    # Exercise functionality
    B_wedge_C = B_vector ∧ C_vector
    @test B_wedge_C ≈ expected_B_wedge_C
    @test B_wedge_C == outer(B, C)
end
