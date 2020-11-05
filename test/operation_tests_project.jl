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
using LinearAlgebra: norm, ⋅
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

# ------ B::Blade

@testset "project(B::Blade, C::Blade)" begin
    # --- Preparations

    test_dim = 15

    # --- Tests

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
    @test norm(projection_coefficients)^2 ≈ grade(B)

    # Check that norm(project(B, C)) TODO

    # ------ grade(B) > grade(C)

    B = Blade(rand(test_dim, 10))
    C = Blade(rand(test_dim, 3))

    @test iszero(project(B, C))

    # ------ dim(B) != dim(C)

    B = Blade(rand(test_dim, 3))
    C = Blade(rand(test_dim + 1, 4))
    @test_throws DimensionMismatch project(B, C)

    # ------ Check consistency with project(v::Vector, B::Blade)

    v = Vector(rand(test_dim))
    B = Blade(rand(test_dim, 3))
    @test project(v, B) ≈ project(Blade(v), B)
end

@testset "project(B::Blade, C::Pseudoscalar)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = 5
    test_vectors = rand(test_dim, 3)

    # --- Tests

    # dim(B) == dim(C)
    B = Blade(test_vectors)
    C = Pseudoscalar(test_dim, test_value)
    @test project(B, C) === B

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value)
    C = Blade(test_vectors)
    @test_throws DimensionMismatch project(B, C)
end

@testset "project(B::Blade, C::Scalar)" begin
    test_dim = 5
    B = Blade(rand(test_dim, 3))

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Blade, C::Vector)" begin
    # --- Preparations

    test_dim = 15
    C = rand(test_dim)

    # --- Tests

    # grade(B) > 1, return_blade == true
    B = Blade(rand(test_dim, 5))
    projection_vectors = basis(B) * transpose(basis(B)) * C
    @test iszero(project(B, C))

    # grade(B) > 1, return_blade == false
    @test project(B, C, return_blade=false) == 0

    # grade(B) == 1, return_blade == true
    B = Blade(rand(test_dim, 1))
    projection_vectors = (C ⋅ basis(B)) * C
    @test project(B, C) ≈ Blade(projection_vectors)

    # grade(B) == 1, return_blade == false
    @test project(B, C, return_blade=false) ≈ projection_vectors

    # dim(B) != length(C)
    B = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch project(B, C)
end

# ------ B::Pseudoscalar

@testset "project(B::Pseudoscalar, C::Blade)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = 10
    test_vectors = rand(test_dim, 3)

    # --- Tests

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value)
    C = Blade(test_vectors)
    @test iszero(project(B, C))

    # dim(B) != dim(C)
    B = Blade(test_vectors)
    C = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch project(B, C)
end

@testset "project(B::Pseudoscalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    for test_dim in 5:8
        # dim(B) == dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_proj_C = project(B, C)
        @test B_proj_C == B

        # dim(B) != dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim + 1, test_value_2)
        @test_throws DimensionMismatch project(B, C)
    end
end

@testset "project(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Pseudoscalar, C::Vector)" begin
    # --- Preparations

    test_dim = 10

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = Vector(rand(test_dim))

    # --- Tests

    # dim(B) == length(C), return_blade == true
    @test iszero(project(B, C))

    # dim(B) == length(C), return_blade == false
    @test project(B, C, return_blade=false) == 0

    # dim(B) != length(C)
    B = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch project(B, C)
end

# ------ B::Scalar

@testset "project(B::Scalar, C::Blade)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    test_dim = 5
    C = Blade(rand(test_dim, 3))

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value
end

@testset "project(B::Scalar, C::Pseudoscalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_dim = 10
    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value_1

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value_1
end

@testset "project(B::Scalar, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value_1

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value_1
end

@testset "project(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    # return_blade == true
    @test project(B, C) === B

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value
end

@testset "project(B::Scalar, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Zero()

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Scalar, C::Real)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value_1

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value_1
end


@testset "project(B::Scalar, C::Vector)" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    test_dim = 15
    C = rand(test_dim)

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value
end

# ------ B::One

@testset "project(B::One, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    # return_blade == true
    @test project(B, C) === B

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 1
end

@testset "project(B::One, C::One)" begin
    B = One()
    C = One()

    # return_blade == true
    @test isone(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 1
end

@testset "project(B::One, C::Zero)" begin
    B = One()
    C = Zero()

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::One, C::Real)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    # return_blade == treu
    @test project(B, C) === B

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 1
end

@testset "project(B::One, C::Vector)" begin
    # return_blade == true
    B = One()
    C = Vector(rand(10))

    # return_blade == true
    @test project(B, C) === B

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 1
end

# ------ B::Zero

@testset "project(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Zero, C::One)" begin
    B = Zero()
    C = One()

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Zero, C::Real)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Zero, C::Vector)" begin
    B = Zero()
    C = Vector(rand(10))

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

# ------ B::Real

@testset "project(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    # return_blade == true
    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value_1

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == test_value_1
end

@testset "project(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_proj_C = project(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == test_value
end

@testset "project(B::Real, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    @test iszero(project(B, C))
end

# ------ B::Vector

@testset "project(B::Vector, C::Blade)" begin
    # --- Preparations

    test_dim = 15
    B = rand(test_dim)

    # --- Tests

    # grade(C) > 1, return_blade == true
    C = Blade(rand(test_dim, 5))
    projection_vectors = basis(C) * transpose(basis(C)) * B
    @test project(B, C) ≈ Blade(projection_vectors)

    # grade(C) > 1, return_blade == false
    @test project(B, C, return_blade=false) ≈ projection_vectors

    # grade(C) == 1, return_blade == true
    C = Blade(rand(test_dim, 1))
    projection_vectors = B ⋅ basis(C) * basis(C)
    @test project(B, C) ≈ Blade(projection_vectors)

    # grade(C) == 1, return_blade == false
    projection_vectors = (B ⋅ basis(C)) * basis(C)
    @test project(B, C, return_blade=false) ≈ projection_vectors

    # length(B) != dim(C)
    C = Blade(rand(test_dim + 1, 3))
    @test_throws DimensionMismatch project(B, C)
end

@testset "project(B::Vector, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10
    B = Vector(rand(test_dim))

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Pseudoscalar(test_dim, test_value)

    # --- Tests

    # length(B) == dim(C), return_blade == true
    @test project(B, C) == Blade(B)

    # length(B) == dim(C), return_blade == false
    @test project(B, C, return_blade=false) == B

    # length(B) != dim(C)
    B = Vector(rand(test_dim + 1))
    @test_throws DimensionMismatch project(B, C)
end

@testset "project(B::Vector, C::Scalar)" begin
    # --- Preparations

    test_dim = 15
    B = rand(test_dim)

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    # return_blade == true
    @test iszero(project(B, C))

    # return_blade == false
    B_proj_C = project(B, C, return_blade=false)
    @test B_proj_C isa Real
    @test B_proj_C == 0
end

@testset "project(B::Vector, C::One)" begin
    # return_blade == true
    B = Vector(rand(10))
    C = One()
    @test iszero(project(B, C))

    # return_blade == false
    @test project(B, C, return_blade=false) == 0
end

@testset "project(B::Vector, C::Zero)" begin
    B = Vector(rand(10))
    C = Zero()
    @test iszero(project(B, C))
end
