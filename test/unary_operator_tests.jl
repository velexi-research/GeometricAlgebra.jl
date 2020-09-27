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


# --- -(x)

@testset "-(x) tests: x::Blade" begin
    # Preparations
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)

    # x::Blade
    negative_B = Blade(B, volume=-volume(B))
    @test -B == negative_B

    @test -negative_B == B
end

@testset "-(x) tests: x::Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Scalar{<:AbstractFloat}
    for value_type in subtypes(AbstractFloat)
        B = Scalar(value_type(value))
        expected_result = Scalar(-value_type(value))
        @test -B == expected_result
    end

    # x::Scalar{<:Signed}
    int_value::Int = 3
    for value_type in subtypes(Signed)
        B = Scalar(value_type(int_value))
        expected_result = Scalar(-value_type(int_value))
        @test -B == expected_result
    end
end

@testset "-(x) tests: x::Pseudoscalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value
    dim = 10

    # x::Pseudoscalar{<:AbstractFloat}
    for value_type in subtypes(AbstractFloat)
        B = Pseudoscalar(dim, value_type(value))
        expected_result = Pseudoscalar(dim, -value_type(value))
        @test -B == expected_result
    end

    # x::Pseudoscalar{<:Signed}
    int_value::Int = 3
    for value_type in subtypes(Signed)
        B = Pseudoscalar(dim, value_type(int_value))
        expected_result = Pseudoscalar(dim, -value_type(int_value))
        @test -B == expected_result
    end
end

# --- reciprocal(x)

@testset "reciprocal(x) tests: x::Blade" begin
    for precision_type in subtypes(AbstractFloat)
        # mod(grade, 4) == 1
        vectors = Vector{precision_type}([3; 4; 0; 0; 0])
        B = Blade(vectors)
        expected_reciprocal = Blade(B, volume=1 / precision_type(5))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(grade, 4) == 2
        vectors = Matrix{precision_type}([3 3; 4 4; 0 1; 0 0; 0 0])
        B = Blade(vectors)
        expected_reciprocal = Blade(B, volume=-1 / precision_type(5))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(grade, 4) == 3
        vectors = Matrix{precision_type}([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
        B = Blade(vectors)
        expected_reciprocal = Blade(B, volume=-1 / precision_type(5))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(grade, 4) == 0
        vectors = Matrix{precision_type}(
            [3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
        B = Blade(vectors)
        expected_reciprocal = Blade(B, volume=1 / precision_type(5))
        @test reciprocal(B) ≈ expected_reciprocal
    end
end

@testset "reciprocal(x) tests: x::Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Scalar
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_value = precision_type(value)

        # value > 0
        S = Scalar{precision_type}(converted_value)
        @test reciprocal(S) == Scalar{precision_type}(1 / converted_value)

        # value < 0
        negative_value = -(abs(converted_value))
        S = Scalar{precision_type}(negative_value)
        @test reciprocal(S) == Scalar{precision_type}(1 / negative_value)

        # value = 0
        S = zero(Scalar{precision_type})
        @test reciprocal(S) == Scalar{precision_type}(Inf)

        # value == 1
        S = one(Scalar{precision_type})
        @test reciprocal(S) == S

        # value = Inf
        S = Scalar{precision_type}(Inf)
        @test reciprocal(S) == zero(Scalar{precision_type})

        # value = -Inf
        S = Scalar{precision_type}(-Inf)
        @test reciprocal(S) == zero(Scalar{precision_type})
    end
end

@testset "reciprocal(x) tests: x::Pseudoscalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Pseudoscalar
    for precision_type in subtypes(AbstractFloat)
        # mod(dim, 4) == 1
        dim = 5
        B = Pseudoscalar(dim, value)
        expected_reciprocal = Pseudoscalar(dim, 1 / precision_type(value))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(dim, 4) == 2
        dim = 6
        B = Pseudoscalar(dim, value)
        expected_reciprocal = Pseudoscalar(dim, -1 / precision_type(value))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(dim, 4) == 3
        dim = 7
        B = Pseudoscalar(dim, value)
        expected_reciprocal = Pseudoscalar(dim, -1 / precision_type(value))
        @test reciprocal(B) ≈ expected_reciprocal

        # mod(dim, 4) == 0
        dim = 8
        B = Pseudoscalar(dim, value)
        expected_reciprocal = Pseudoscalar(dim, 1 / precision_type(value))
        @test reciprocal(B) ≈ expected_reciprocal
    end
end
