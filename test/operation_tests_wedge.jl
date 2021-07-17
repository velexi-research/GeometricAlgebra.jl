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
using LinearAlgebra: norm
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

#=
# --- Tests

# ------ M::Multivector

@testset "wedge(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "wedge(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "wedge(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Blade, C::Blade)" begin
    # --- Preparations

    B_vectors = hcat([1; 1; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vectors = hcat([0; 1; 3; 0; 0],
                     [0; 0; 0; 4; 0])
    C = Blade(C_vectors)

    # --- Tests

    # ------ dim(B) == dim(C)

    # grade(B) + grade(C) <= dim(B)
    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(hcat(B_vectors, C_vectors))
    @test B ∧ C == B_wedge_C

    # grade(B) + grade(C) > dim(B)
    C_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 1; 0; 0; 0],
                     [0; 0; 1; 0; 0],
                     [0; 0; 0; 1; 0])
    C = Blade(C_vectors)
    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)

    # ------ dim(B) != dim(C)

    C = Blade(rand(dim(B) + 1, 2))
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Blade, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10

    test_vectors = rand(test_dim, 3)

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Test

    # dim(B) == dim(C)
    B = Blade(test_vectors)
    C = Pseudoscalar(test_dim, test_value)

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)

    # dim(B) != dim(C)
    B = Blade(test_vectors)
    C = Pseudoscalar(test_dim + 1, test_value)
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Blade, C::Scalar)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(basis(B), volume=value(C) * volume(B))
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Blade, C::One)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    C = One()

    B_wedge_C = wedge(B, C)
    @test B_wedge_C === B
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Blade, C::Zero)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    C = Zero()

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Blade, C::Real)" begin
    test_dim = 10
    B = Blade(rand(test_dim, 3))

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(basis(B), volume=C * volume(B))
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Blade, C::Vector)" begin
    # --- Preparations

    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C = [0; 1; 3; 0; 0]

    # --- Tests

    # dim(B) == dim(C)
    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(hcat(B_vectors, C))
    @test B ∧ C == B_wedge_C

    # dim(B) != dim(C)
    C = Vector(rand(dim(B) + 1))
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

# ------ B::Pseudoscalar

@testset "wedge(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Pseudoscalar, C::Blade)" begin
    # --- Preparations

    test_dim = 10

    test_vectors = rand(test_dim, 3)

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Test

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value)
    C = Blade(test_vectors)

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value)
    C = Blade(test_vectors)
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Pseudoscalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Pseudoscalar, C::Scalar)" begin
    # --- Preparations

    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = One()

    B_wedge_C = wedge(B, C)
    @test B_wedge_C === B
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Pseudoscalar(test_dim, test_value)

    C = Zero()

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Pseudoscalar, C::Real)" begin
    # --- Preparations

    test_dim = 10
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Pseudoscalar(test_dim, test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Pseudoscalar, C::Vector)" begin
    # --- Preparations

    test_dim = 10

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # dim(B) == dim(C)
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value)
        C = Vector(rand(test_dim))

        @test iszero(wedge(B, C))
        @test iszero(B ∧ C)
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim + 1, test_value)
    C = Vector(rand(test_dim))

    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

# ------ B::Scalar

@testset "wedge(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Scalar, C::Blade)" begin
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(basis(C), volume=value(B) * volume(C))
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Scalar, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_dim = 10
    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Scalar, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test B_wedge_C == test_value_1 * test_value_2
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Scalar, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = One()

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test B_wedge_C == test_value
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Scalar, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Zero()

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Scalar, C::Real)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = Scalar(test_value_1)

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = test_value_2

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test B_wedge_C == test_value_1 * test_value_2
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Scalar, C::Vector)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = Scalar(test_value)

    C = Vector(rand(5))

    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(C, volume=norm(C) * value(B))
    @test B ∧ C == B_wedge_C
end

# ------ B::One

@testset "wedge(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::One::, C::Blade)" begin
    B = One()

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    B_wedge_C = wedge(B, C)
    @test B_wedge_C === C
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::One::, C::Pseudoscalar)" begin
    B = One()

    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Pseudoscalar(test_dim, test_value)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C === C
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::One::, C::Scalar)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    @test wedge(B, C) isa AbstractScalar
    @test wedge(B, C) == test_value
    @test B ∧ C isa AbstractScalar
    @test B ∧ C == test_value
end

@testset "wedge(B::One, C::One)" begin
    B = One()
    C = One()
    @test isone(wedge(B, C))
    @test isone(B ∧ C)
end

@testset "wedge(B::One, C::Zero)" begin
    B = One()
    C = Zero()
    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::One, C::Real)" begin
    B = One()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test wedge(B, C) == C
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::One, C::Vector)" begin
    B = One()
    C = Vector(rand(5))

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Blade(C)
    @test B ∧ C == B_wedge_C
end

# ------ B::Zero

@testset "wedge(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Zero, C::Blade)" begin
    B = Zero()

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::Pseudoscalar)" begin
    B = Zero()

    test_dim = 10
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Pseudoscalar(test_dim, test_value)

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::Scalar)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::One)" begin
    B = Zero()
    C = One()
    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::Zero)" begin
    B = Zero()
    C = Zero()
    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::Real)" begin
    B = Zero()

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = test_value

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Zero, C::Vector)" begin
    B = Zero()
    C = Vector(rand(5))

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

# ------ B::Real

@testset "wedge(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Real, C::Blade)" begin
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    test_dim = 10
    C = Blade(rand(test_dim, 3))

    # Exercise functionality and check results
    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(basis(C), volume=B * volume(C))
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Real, C::Pseudoscalar)" begin
    # --- Preparations

    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_dim = 10
    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Pseudoscalar(test_dim, test_value_2)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Real, C::Scalar)" begin
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1
    B = test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2
    C = Scalar(test_value_2)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test B_wedge_C == test_value_1 * test_value_2
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Real, C::One)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = One()

    B_wedge_C = wedge(B, C)
    @test B_wedge_C isa AbstractScalar
    @test B_wedge_C == test_value
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Real, C::Zero)" begin
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    B = test_value

    C = Zero()

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

# ------ B::Vector

@testset "wedge(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "wedge(B::Vector, C::Blade)" begin
    # --- Preparations

    B = [0; 1; 3; 0; 0]

    C_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    C = Blade(C_vectors)

    # --- Tests

    # dim(B) == dim(C)
    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(hcat(B, C_vectors))
    @test B ∧ C == B_wedge_C

    # dim(B) != dim(C)
    B = Vector(rand(dim(C) + 1))
    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Vector, C::Pseudoscalar)" begin
    # --- Preparations

    test_dim = 10

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # dim(B) == dim(C)
    for test_dim in 5:8
        B = Vector(rand(test_dim))
        C = Pseudoscalar(test_dim, test_value)

        @test iszero(wedge(B, C))
        @test iszero(B ∧ C)
    end

    # dim(B) != dim(C)
    B = Vector(rand(test_dim))
    C = Pseudoscalar(test_dim + 1, test_value)

    @test_throws DimensionMismatch wedge(B, C)
    @test_throws DimensionMismatch B ∧ C
end

@testset "wedge(B::Vector, C::Scalar)" begin
    B = Vector(rand(5))

    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value
    C = Scalar(test_value)

    B_wedge_C = wedge(B, C)
    @test B_wedge_C ≈ Blade(B, volume=norm(B) * value(C))
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Vector, C::One)" begin
    B = Vector(rand(5))
    C = One()

    B_wedge_C = wedge(B, C)
    @test B_wedge_C == Blade(B)
    @test B ∧ C == B_wedge_C
end

@testset "wedge(B::Vector, C::Zero)" begin
    B = Vector(rand(5))
    C = Zero()

    @test iszero(wedge(B, C))
    @test iszero(B ∧ C)
end

@testset "wedge(B::Vector, C::Vector)" begin
    # --- Preparations

    B_vector = Vector{Float64}([1, 2, 3, 0, 0])
    B = Blade(B_vector)

    C_vector = [0; 1; 3; 0; 0]
    C = Blade(C_vector)

    # --- Tests

    # dim(B) == dim(C)
    B_wedge_C = wedge(B_vector, C_vector)
    @test B_wedge_C ≈ Blade(hcat(B_vector, C_vector))
    @test B_vector ∧ C_vector == B_wedge_C

    # dim(B) != dim(C)
    C_vector = Vector(rand(dim(B) + 1))
    @test_throws DimensionMismatch wedge(B_vector, C_vector)
    @test_throws DimensionMismatch B_vector ∧ C_vector
end
=#
