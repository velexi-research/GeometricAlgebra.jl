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
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
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
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2
end
