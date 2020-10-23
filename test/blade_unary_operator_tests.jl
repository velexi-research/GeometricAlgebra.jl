"""
Unit tests for AbstractBlade unary operators.

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


# --- -(B)

@testset "-(B): B::Blade" begin
    # Preparations
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)

    # B::Blade
    negative_B = Blade(B, volume=-volume(B))
    @test -B == negative_B

    @test -negative_B == B
end

@testset "-(B): B::Scalar" begin
    # Preparations
    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # B::Scalar{<:AbstractFloat}
    for value_type in subtypes(AbstractFloat)
        B = Scalar(value_type(test_value))
        expected_result = Scalar(-value_type(test_value))
        @test -B == expected_result
    end

    # B::Scalar{<:Signed}
    int_value::Int = 3
    for value_type in subtypes(Signed)
        B = Scalar(value_type(int_value))
        expected_result = Scalar(-value_type(int_value))
        @test -B == expected_result
    end
end

@testset "-(B): B::Pseudoscalar" begin
    # Preparations
    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value
    dim = 10

    # B::Pseudoscalar{<:AbstractFloat}
    for value_type in subtypes(AbstractFloat)
        B = Pseudoscalar(dim, value_type(test_value))
        expected_result = Pseudoscalar(dim, -value_type(test_value))
        @test -B == expected_result
    end

    # B::Pseudoscalar{<:Signed}
    int_value::Int = 3
    for value_type in subtypes(Signed)
        B = Pseudoscalar(dim, value_type(int_value))
        expected_result = Pseudoscalar(dim, -value_type(int_value))
        @test -B == expected_result
    end
end

# --- reciprocal(B)

@testset "reciprocal(B): B::Blade" begin
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
    @test_skip B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test_skip B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test_skip B * reciprocal(B) ≈ 1
end


@testset "reciprocal(B): B::Scalar" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    # value > 0
    B = Scalar(test_value)
    @test reciprocal(B) == Scalar(1 / test_value)
    @test B * reciprocal(B) ≈ 1

    # value < 0
    negative_value = -(abs(test_value))
    B = Scalar(negative_value)
    @test reciprocal(B) == Scalar(1 / negative_value)
    @test B * reciprocal(B) ≈ 1

    # value = 0
    B = zero(B)
    @test reciprocal(B) == Scalar(Inf)

    # value == 1
    B = one(B)
    @test reciprocal(B) == B
    @test B * reciprocal(B) ≈ 1

    # value = Inf
    B = Scalar(Inf)
    @test reciprocal(B) == zero(B)

    # value = -Inf
    B = Scalar(-Inf)
    @test reciprocal(B) == zero(B)
end


@testset "reciprocal(B): B::Pseudoscalar" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    # mod(dim, 4) == 1
    dim = 5
    B = Pseudoscalar(dim, test_value)
    expected_reciprocal = Pseudoscalar(dim, 1 / test_value)
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(dim, 4) == 2
    dim = 6
    B = Pseudoscalar(dim, test_value)
    expected_reciprocal = Pseudoscalar(dim, -1 / test_value)
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(dim, 4) == 3
    dim = 7
    B = Pseudoscalar(dim, test_value)
    expected_reciprocal = Pseudoscalar(dim, -1 / test_value)
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(dim, 4) == 0
    dim = 8
    B = Pseudoscalar(dim, test_value)
    expected_reciprocal = Pseudoscalar(dim, 1 / test_value)
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1
end


# --- reverse(B)

@testset "reverse(B): B::Blade" begin
    # mod(grade, 4) == 1
    vectors = Vector([3; 4; 0; 0; 0])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 2
    vectors = Matrix([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test_skip B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test_skip B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    @test reverse(B) === B
    @test_skip B * reverse(B) ≈ norm(B)^2
end


@testset "reverse(B): B::Scalar" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    # value != 0
    B = Scalar(test_value)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # value = 0
    B = zero(B)
    @test reverse(B) === B

    # value == 1
    B = one(B)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # value = Inf
    B = Scalar(Inf)
    @test reverse(B) === B

    # value = -Inf
    B = Scalar(-Inf)
    @test reverse(B) === B
end


@testset "reverse(B): B::Pseudoscalar" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # mod(dim, 4) == 1
    dim = 5
    B = Pseudoscalar(dim, test_value)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(dim, 4) == 2
    dim = 6
    B = Pseudoscalar(dim, test_value)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(dim, 4) == 3
    dim = 7
    B = Pseudoscalar(dim, test_value)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(dim, 4) == 0
    dim = 8
    B = Pseudoscalar(dim, test_value)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2
end
