"""
Unit tests for Blade operators.

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
import LinearAlgebra
using Test

# GeometricAlgebra.jl
using GeometricAlgebra


# --- ==(x, y)

@testset "==(x, y): Blade" begin
    # --- Preparations

    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # --- Test functionality

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors))
            @test B1 ≈ B2
        end
    end
end

@testset "==(x, y): Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    float64_or_bigfloat = (Float64, BigFloat)

    # S1::Scalar, S2::Scalar
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        S = Scalar(converted_value)
        @test S == converted_value
        @test converted_value == S
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            S1 = Scalar(precision_type1(value))
            S2 = Scalar(precision_type2(value))
            if precision_type1 == precision_type2
                @test S1 == S2
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test S1 == S2
            else
                @test S1 != S2
            end
        end
    end

    # S1::Scalar, S2::Real
    # S1::Real, S2::Scalar
    for precision_type in subtypes(AbstractFloat)
        # S1::Scalar, S2::AbstractFloat
        # S1::AbstractFloat, S2::Scalar
        for value_type in subtypes(AbstractFloat)
            if precision_type == value_type
                @test Scalar(precision_type(value)) == value_type(value)
                @test value_type(value) == Scalar(precision_type(value))
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test Scalar(precision_type(value)) == value_type(value)
                @test value_type(value) == Scalar(precision_type(value))
            else
                @test Scalar(precision_type(value)) != value_type(value)
                @test value_type(value) != Scalar(precision_type(value))
            end
        end

        # S1::Scalar, S2::Integer
        # S1::Integer, S2::Scalar
        int_value::Int = 3
        for value_type in subtypes(Signed)
            @test Scalar(precision_type(int_value)) == value_type(int_value)
            @test value_type(int_value) == Scalar(precision_type(int_value))
        end

        for value_type in subtypes(Unsigned)
            @test Scalar(precision_type(int_value)) == value_type(int_value)
            @test value_type(int_value) == Scalar(precision_type(int_value))
        end

        @test Scalar(precision_type(true)) == true
        @test true == Scalar(precision_type(true))
        @test Scalar(precision_type(false)) == false
        @test false == Scalar(precision_type(false))
    end

    # S1::Scalar, S2::Zero
    # S1::Zero, S2::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            @test Scalar(precision_type1(0)) == Zero(precision_type2)
            @test Zero(precision_type2) == Scalar(precision_type1(0))
        end
    end

    # S1::Scalar, S2::One
    # S1::One, S2::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            @test Scalar(precision_type1(1)) == One(precision_type2)
            @test One(precision_type2) == Scalar(precision_type1(1))
        end
    end

    # S1::Zero, S2::Zero
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            @test Zero(precision_type1) == Zero(precision_type2)
        end
    end

    # S1::Zero, S2::Real
    # S1::Real, S2::Zero
    for precision_type in subtypes(AbstractFloat)
        @test Zero(precision_type) == 0
        @test 0 == Zero(precision_type)
    end

    # S1::One, S2::One
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            @test One(precision_type1) == One(precision_type2)
        end
    end

    # S1::One, S2::Real
    # S1::Real, S2::One
    for precision_type in subtypes(AbstractFloat)
        @test One(precision_type) == 1
        @test 1 == One(precision_type)
    end
end

# --- ≈(x, y)

@testset "≈(x, y) Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        S = Scalar(converted_value)
        @test S ≈ converted_value
        @test converted_value ≈ S
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            S1 = Scalar(precision_type1(value))
            S2 = Scalar(precision_type2(value))
            @test S1 ≈ S2
        end
    end
end

@testset "-() tests: Blade" begin
    # mod(grade, 4) == 0
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_result = Blade(B, sign=-1)
    @test -B == expected_result
end

@testset "inverse() tests: Blade" begin
    # mod(grade, 4) == 1
    vectors = Vector{BigFloat}([3; 4; 0; 0; 0])
    B = Blade(vectors)
    expected_inverse_norm = 1 / BigFloat(25)
    expected_inverse = Blade(
        [expected_inverse_norm; expected_inverse_norm; 1; 1; 1] .* vectors)
    @test inverse(B) ≈ expected_inverse

    # mod(grade, 4) == 2
    vectors = Matrix{Float64}([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    expected_inverse_norm = 1 / Float64(25)
    expected_inverse = -Blade(
        [expected_inverse_norm; expected_inverse_norm; 1; 1; 1] .* vectors)
    @test inverse(B) ≈ expected_inverse

    # mod(grade, 4) == 3
    vectors = Matrix{Float32}([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    expected_inverse_norm = 1 / Float32(25)
    expected_inverse = -Blade(
        [expected_inverse_norm; expected_inverse_norm; 1; 1; 1] .* vectors)
    @test inverse(B) ≈ expected_inverse

    # mod(grade, 4) == 0
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_inverse_norm = 1 / Float16(25)
    expected_inverse = Blade(
        [expected_inverse_norm; expected_inverse_norm; 1; 1; 1] .* vectors)
    @test inverse(B) ≈ expected_inverse
end
