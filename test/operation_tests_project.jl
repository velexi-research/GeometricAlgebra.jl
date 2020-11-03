"""
Unit tests for the project(x, y) function

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

@testset "project(B, C): B or C isa Blade" begin
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

    projection = project(B, C)
    @test projection isa Blade
    @test dim(projection) == dim(B)

    # Check norm
    norm_projection = norm(B) *
        norm(Blade(basis(C) * transpose(basis(C)) * basis(B)))
    @test norm(projection) ≈ norm_projection

    # Check that project(B, C) is contained in C
    projection_coefficients = transpose(basis(C)) * basis(projection)
    @test LinearAlgebra.norm(projection_coefficients)^2 ≈ grade(B)

    # Check that norm(project(B, C)) TODO

    # ------ grade(B) > grade(C)

    B = Blade(rand(test_dim, 10))
    C = Blade(rand(test_dim, 3))

    @test project(B, C) == zero(B)

    # ------ dim(B) != dim(C)

    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch project(B, C)

    # ------ Check consistency with project(v::Vector, B::Blade)

    v = Vector(rand(test_dim))
    B = Blade(rand(test_dim, 3))
    @test project(v, B) ≈ project(Blade(v), B)

    # --- B::Pseudoscalar, C::Blade
    #     B::Blade, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Blade(test_vectors)

    expected_result = zero(B)
    @test project(B, C) == expected_result

    expected_result = C
    @test project(C, B) == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value_1)
    C = Blade(test_vectors)
    @test_throws DimensionMismatch project(B, C)
    @test_throws DimensionMismatch project(C, B)

    # --- B::Scalar, C::Blade
    #     B::Blade, C::Scalar

    # Preparations
    B = Scalar(test_value_1)
    C = Blade(test_vectors)

    # Exercise functionality and check results
    expected_result = Scalar(value(B))
    @test project(B, C) == expected_result

    expected_result = zero(C)
    @test project(C, B) == expected_result

    # ------ project(v::Vector{<:Real}, B::Scalar)

    # return_blade == true
    B = Scalar(test_value)
    @test project(v, B) == zero(B)

    # return_blade == false
    @test project(v, B, return_blade=false) == 0

    # ------ project(B::Scalar, v::Vector{<:Real})

    # return_blade == true
    B = Scalar(test_value)
    @test project(B, v) === B

    # return_blade == false
    B = Scalar(test_value)
    @test project(B, v, return_blade=false) == value(B)

    # ------ project(v::Vector{<:Real}, B::Pseudoscalar)
    #        project(B::Pseudoscalar, v::Vector{<:Real})

    # length(v) == dim(B), return_blade == true
    B = Pseudoscalar(test_dim, test_value)
    @test project(v, B) == Blade(v)
    @test project(B, v) == zero(B)

    # length(v) == dim(B), return_blade == false
    @test project(v, B, return_blade=false) == v
    @test project(B, v, return_blade=false) == 0

    # length(v) != dim(B)
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch project(v, B)
    @test_throws DimensionMismatch project(B, v)

    # ------ project(v::Vector{<:Real}, B::Blade)
    #        project(B::Blade, v::Vector{<:Real})

    # grade(B) > 1, return_blade == true
    B = Blade(rand(test_dim, 5))
    projection_vectors = basis(B) * transpose(basis(B)) * v
    @test project(v, B) ≈ Blade(projection_vectors)
    @test project(B, v) == zero(B)

    # grade(B) > 1, return_blade == false
    @test project(v, B, return_blade=false) ≈ projection_vectors
    @test project(B, v, return_blade=false) == 0

    # grade(B) == 1, return_blade == true
    B = Blade(rand(test_dim, 1))
    projection_vectors = LinearAlgebra.dot(v, basis(B)) * basis(B)
    @test project(v, B) ≈ Blade(projection_vectors)

    projection_vectors = LinearAlgebra.dot(v, basis(B)) * v
    @test project(B, v) ≈ Blade(projection_vectors)

    # grade(B) == 1, return_blade == false
    projection_vectors = LinearAlgebra.dot(v, basis(B)) * basis(B)
    @test project(v, B, return_blade=false) ≈ projection_vectors

    projection_vectors = LinearAlgebra.dot(v, basis(B)) * v
    @test project(B, v, return_blade=false) ≈ projection_vectors

    # length(v) != dim(B)
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch project(v, B)
    @test_throws DimensionMismatch project(B, v)
end

@testset "project(B, C): B or C isa Pseudoscalar" begin
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

    for test_dim in 5:8
        # dim(B) == dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_proj_C = project(B, C)
        expected_result = B
        @test B_proj_C == expected_result

        # dim(B) != dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim + 1, test_value_2)
        @test_throws DimensionMismatch project(B, C)
    end

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = test_value_1
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == expected_result

    expected_result = zero(C)
    C_proj_B = project(C, B)
    @test C_proj_B === expected_result
end

@testset "project(B, C): B or C isa Scalar" begin
    # --- Preparations

    # Test dimension
    test_dim = 15

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # Test vectors
    test_vector = rand(test_dim)

    # --- B::Scalar, C::Scalar
    B = Scalar(test_value_1)
    C = Scalar(test_value_2)

    B_proj_C = project(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    # --- B::Scalar, C::One
    #     B::One, C::Scalar
    B = Scalar(test_value_1)
    C = One()

    B_proj_C = project(B, C)
    expected_result = B
    @test B_proj_C == expected_result

    C_proj_B = project(C, B)
    expected_result = C
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Zero
    #     B::Zero, C::Scalar
    B = Scalar(test_value_1)
    C = Zero()

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # --- B::Scalar, C::Real
    #     B::Real, C::Scalar
    B = Scalar(test_value_1)
    C = test_value_2

    B_proj_C = project(B, C)
    expected_result = test_value_1
    @test B_proj_C isa Scalar
    @test B_proj_C == expected_result

    C_proj_B = project(C, B)
    expected_result = test_value_2
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result

    # --- B::Scalar, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::Scalar

    # return_blade == true
    B = Scalar(test_value_1)
    C = test_vector
    @test project(B, C) == B
    @test project(C, B) === zero(B)

    # return_blade == false
    @test project(C, B, return_blade=false) == 0
    @test project(B, C, return_blade=false) == value(B)
end

@testset "project(B, C): B or C isa One" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_proj_C = project(B, C)
    expected_result = One()
    @test B_proj_C === expected_result

    # B::One, C::Zero
    # B::Zero, C::One
    B = One()
    C = Zero()

    expected_result = Zero()

    B_proj_C = project(B, C)
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    @test C_proj_B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_proj_C = project(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = C
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result

    # --- B::One, C::Vector{<:Real}
    #     B::Vector{<:Real}, C::One

    # return_blade == true
    B = One()
    C = Vector(rand(10))
    @test project(B, C) === B
    @test project(C, B) == zero(B)

    # return_blade == false
    @test project(C, B, return_blade=false) == 0
    @test project(B, C, return_blade=false) == value(B)
end

@testset "project(B, C): B or C isa Zero" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::Zero, C::Zero
    B = Zero()
    C = Zero()

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    # B::Zero, C::Real
    # B::Real, C::Zero
    B = Zero()
    C = test_value

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result

    # B::Zero, C::Vector
    # B::Vector, C::Zero
    B = Zero()
    C = Vector(rand(10))

    B_proj_C = project(B, C)
    expected_result = Zero()
    @test B_proj_C === expected_result

    C_proj_B = project(C, B)
    expected_result = Zero()
    @test C_proj_B === expected_result
end
