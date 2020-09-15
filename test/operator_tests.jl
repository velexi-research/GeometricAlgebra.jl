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
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # x::Blade, y::Blade
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors))
            if precision_type1 == precision_type2
                @test B1 == B2
            else
                @test B1 != B2
            end
        end
    end
end

@testset "==(x, y), !=(x, y): Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    float64_or_bigfloat = (Float64, BigFloat)

    # x::Scalar, y::Scalar
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        B = Scalar(converted_value)

        # ==(x, y)
        @test B == converted_value
        @test converted_value == B

        # !=(x,y)
        @test B != 2 * converted_value
        @test 2 * converted_value != B
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(x, y)
            B1 = Scalar(precision_type1(value))
            B2 = Scalar(precision_type2(value))
            if precision_type1 == precision_type2
                @test B1 == B2
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B1 == B2
            else
                @test B1 != B2
            end

            # !=(x,y)
            B1 = Scalar(precision_type1(value))
            B2 = Scalar(precision_type2(2 * value))
            @test B1 != B2
        end
    end

    # x::Scalar, y::Real
    # x::Real, y::Scalar
    for precision_type in subtypes(AbstractFloat)
        # x::Scalar, y::AbstractFloat
        # x::AbstractFloat, y::Scalar
        for value_type in subtypes(AbstractFloat)
            # ==(x, y)
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

            # !=(x,y)
            @test Scalar(precision_type(value)) != value_type(2 * value)
            @test value_type(2 * value) != Scalar(precision_type(value))
        end

        # x::Scalar, y::Integer
        # x::Integer, y::Scalar
        int_value::Int = 3
        for value_type in subtypes(Signed)
            # ==(x, y)
            @test Scalar(precision_type(int_value)) == value_type(int_value)
            @test value_type(int_value) == Scalar(precision_type(int_value))

            # !=(x,y)
            @test Scalar(precision_type(int_value)) != value_type(2 * int_value)
            @test value_type(2 * int_value) != Scalar(precision_type(int_value))
        end

        for value_type in subtypes(Unsigned)
            # ==(x, y)
            @test Scalar(precision_type(int_value)) == value_type(int_value)
            @test value_type(int_value) == Scalar(precision_type(int_value))

            # !=(x,y)
            @test Scalar(precision_type(int_value)) != value_type(2 * int_value)
            @test value_type(2 * int_value) != Scalar(precision_type(int_value))
        end

        # Bool
        @test Scalar(precision_type(true)) == 1
        @test 1 == Scalar(precision_type(true))
        @test Scalar(precision_type(true)) != 0
        @test false != Scalar(precision_type(true))

        @test Scalar(precision_type(false)) == 0
        @test 0 == Scalar(precision_type(false))
        @test Scalar(precision_type(false)) != 1
        @test 1 != Scalar(precision_type(false))
    end

    # x::Scalar, y::Zero
    # x::Zero, y::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # !=(x,y)
            @test Scalar(precision_type1(1)) != Zero(precision_type2)
            @test Zero(precision_type2) != Scalar(precision_type1(1))
        end
    end

    # x::Scalar, y::One
    # x::One, y::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(x, y)
            @test Scalar(precision_type1(1)) == One(precision_type2)
            @test One(precision_type2) == Scalar(precision_type1(1))

            # ==(x, y)
            @test Scalar(precision_type1(5)) != One(precision_type2)
            @test One(precision_type2) != Scalar(precision_type1(5))
        end
    end

    # x::Zero, y::Zero
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(x, y)
            @test Zero(precision_type1) == Zero(precision_type2)
        end
    end

    # x::Zero, y::Real
    # x::Real, y::Zero
    for precision_type in subtypes(AbstractFloat)
        # x::Zero, y::AbstractFloat
        # x::AbstractFloat, y::Zero
        for value_type in subtypes(AbstractFloat)
            # ==(x, y)
            @test Zero(precision_type) == value_type(0)
            @test value_type(0) == Zero(precision_type)

            # !=(x, y)
            @test Zero(precision_type) != value_type(5)
            @test value_type(5) != Zero(precision_type)
        end

        # x::Zero, y::Integer
        # x::Integer, y::Zero
        for value_type in subtypes(Signed)
            # ==(x, y)
            @test Zero(precision_type) == value_type(0)
            @test value_type(0) == Zero(precision_type)

            # !=(x, y)
            @test Zero(precision_type) != value_type(5)
            @test value_type(5) != Zero(precision_type)
        end

        for value_type in subtypes(Unsigned)
            # ==(x, y)
            @test Zero(precision_type) == value_type(0)
            @test value_type(0) == Zero(precision_type)

            # !=(x, y)
            @test Zero(precision_type) != value_type(5)
            @test value_type(5) != Zero(precision_type)
        end

        # Bool
        @test Zero(precision_type) == 0
        @test 0 == Zero(precision_type)
        @test Zero(precision_type) != 1
        @test 1 != Zero(precision_type)
    end

    # x::One, y::One
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(x, y)
            @test One(precision_type1) == One(precision_type2)
        end
    end

    # x::One, y::Real
    # x::Real, y::One
    for precision_type in subtypes(AbstractFloat)
        # x::One, y::AbstractFloat
        # x::AbstractFloat, y::One
        for value_type in subtypes(AbstractFloat)
            # ==(x, y)
            @test One(precision_type) == value_type(1)
            @test value_type(1) == One(precision_type)

            # !=(x, y)
            @test One(precision_type) != value_type(3)
            @test value_type(3) != One(precision_type)
        end

        # x::One, y::Integer
        # x::Integer, y::One
        for value_type in subtypes(Signed)
            # ==(x, y)
            @test One(precision_type) == value_type(1)
            @test value_type(1) == One(precision_type)

            # !=(x, y)
            @test One(precision_type) != value_type(5)
            @test value_type(5) != One(precision_type)
        end

        for value_type in subtypes(Unsigned)
            # ==(x, y)
            @test One(precision_type) == value_type(1)
            @test value_type(1) == One(precision_type)

            # !=(x, y)
            @test One(precision_type) != value_type(5)
            @test value_type(5) != One(precision_type)
        end

        # Bool
        @test One(precision_type) == 1
        @test 1 == One(precision_type)
        @test One(precision_type) != 0
        @test 0 != One(precision_type)
    end

    # x::Zero, y::One
    # x::One, y::Zero
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # !=(x, y)
            @test One(precision_type1) != Zero(precision_type2)
            @test Zero(precision_type1) != One(precision_type2)
        end
    end
end

# --- ≈(x, y)

@testset "≈(x, y): Blade" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # x::Blade, y::Blade
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors))
            @test B1 ≈ B2
        end
    end
end

@testset "≈(x, y): Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Scalar, y::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Scalar(precision_type1(value))
            B2 = Scalar(precision_type2(value))
            @test B1 ≈ B2
        end
    end

    # x::Scalar, y::Real
    # x::Real, y::Scalar
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        B = Scalar(converted_value)
        @test B ≈ converted_value
        @test converted_value ≈ B
    end
end

# --- -(x), opposite(x)

@testset "-(x) tests: Blade" begin
    # Preparations
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)

    # x::Blade
    negative_B = Blade(B, sign=-1)
    @test -B == negative_B
    @test opposite(B) == negative_B

    @test -negative_B == B
    @test opposite(negative_B) == B
end

@testset "-(x) tests: Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Scalar{<:AbstractFloat}
    for value_type in subtypes(AbstractFloat)
        B = Scalar(value_type(value))
        expected_result = Scalar(-value_type(value))
        @test -B == expected_result
        @test opposite(B) == expected_result
    end

    # x::Scalar{<:Signed}
    int_value::Int = 3
    for value_type in subtypes(Signed)
        B = Scalar(value_type(int_value))
        expected_result = Scalar(-value_type(int_value))
        @test -B == expected_result
        @test opposite(B) == expected_result
    end
end

# --- reciprocal(x)

@testset "reciprocal(x) tests: Blade" begin
    # mod(grade, 4) == 1
    vectors = Vector{BigFloat}([3; 4; 0; 0; 0])
    B = Blade(vectors)
    reciprocal_norm = 1 / BigFloat(25)
    expected_reciprocal = Blade(
        [reciprocal_norm; reciprocal_norm; 1; 1; 1] .* vectors)
    @test reciprocal(B) ≈ expected_reciprocal

    # mod(grade, 4) == 2
    vectors = Matrix{Float64}([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    reciprocal_norm = 1 / Float64(25)
    expected_reciprocal = -Blade(
        [reciprocal_norm; reciprocal_norm; 1; 1; 1] .* vectors)
    @test reciprocal(B) ≈ expected_reciprocal

    # mod(grade, 4) == 3
    vectors = Matrix{Float32}([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    reciprocal_norm = 1 / Float32(25)
    expected_reciprocal = -Blade(
        [reciprocal_norm; reciprocal_norm; 1; 1; 1] .* vectors)
    @test reciprocal(B) ≈ expected_reciprocal

    # mod(grade, 4) == 0
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    reciprocal_norm = 1 / Float16(25)
    expected_reciprocal = Blade(
        [reciprocal_norm; reciprocal_norm; 1; 1; 1] .* vectors)
    @test reciprocal(B) ≈ expected_reciprocal
end

@testset "reciprocal(x) tests: Scalar" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    # x::Scalar
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_value = precision_type(value)

        # value > 0
        S = Scalar(converted_value)
        @test reciprocal(S) == Scalar{precision_type}(1 / converted_value)

        # value < 0
        negative_value = -(abs(converted_value))
        S = Scalar(negative_value)
        @test reciprocal(S) == Scalar{precision_type}(1 / negative_value)

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test reciprocal(S) === Zero{precision_type}()

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test reciprocal(S) === Zero{precision_type}()
    end

    # x::Zero
    for precision_type in subtypes(AbstractFloat)
        @test reciprocal(Zero(precision_type)) == Scalar(Inf)
    end

    # x::One
    for precision_type in subtypes(AbstractFloat)
        @test reciprocal(One(precision_type)) === One(precision_type)
    end
end
