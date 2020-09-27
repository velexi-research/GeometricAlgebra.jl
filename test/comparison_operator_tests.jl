"""
Unit tests for AbstractBlade comparison operators.

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


# --- ==(x, y)

@testset "==(x, y): x, y::Blade" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B1) == dim(B2), grade(B1) == grade(B2), volume(B1) == volume(B2)
    # basis(B1) == basis(B2)
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

    # dim(B1) != dim(B2)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors2))
            @test B1 != B2
        end
    end

    # grade(B1) != grade(B2)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors2))
            @test B1 != B2
        end
    end

    # volume(B1) != volume(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors),
                       volume=2*volume(B1))
            @test B1 != B2
        end
    end

    # basis(B1) != basis(B2)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 5; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            B2 = Blade(convert(Array{precision_type2}, vectors2),
                       volume=test_volume)
            @test B1 != B2
        end
    end
end

@testset "==(x, y): x, y::Scalar" begin
    # Preparations
    value = rand()
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
            B = Scalar(precision_type(value))

            # ==(x, y)
            if precision_type == value_type
                @test B == value_type(value)
                @test value_type(value) == B
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test B == value_type(value)
                @test value_type(value) == B
            else
                @test B != value_type(value)
                @test value_type(value) != B
            end

            # !=(x,y)
            @test B != value_type(2 * value)
            @test value_type(2 * value) != B
        end

        # x::Scalar, y::Integer
        # x::Integer, y::Scalar
        int_value::Int = 3
        for value_type in subtypes(Signed)
            B = Scalar(precision_type(int_value))

            # ==(x, y)
            @test B == value_type(int_value)
            @test value_type(int_value) == B

            # !=(x,y)
            @test B != value_type(2 * int_value)
            @test value_type(2 * int_value) != B
        end

        for value_type in subtypes(Unsigned)
            B = Scalar(precision_type(int_value))

            # ==(x, y)
            @test B == value_type(int_value)
            @test value_type(int_value) == B

            # !=(x,y)
            @test B != value_type(2 * int_value)
            @test value_type(2 * int_value) != B
        end

        # Bool
        B = Scalar(precision_type(true))
        @test B == 1
        @test 1 == B
        @test B != 0
        @test 0 != B

        B = Scalar(precision_type(false))
        @test B == 0
        @test 0 == B
        @test B != 1
        @test 1 != B
    end
end

@testset "==(x, y): x, y::Pseudoscalar" begin
    # Preparations
    dim = 10

    value = rand()
    value = rand() > 0.5 ? value : -value

    float64_or_bigfloat = (Float64, BigFloat)

    # dim(B1) == dim(B2), value(B1) == value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim, precision_type2(value))
            if precision_type1 == precision_type2
                @test B1 == B2
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B1 == B2
            else
                @test B1 != B2
            end
        end
    end

    # dim(B1) != dim(B2), value(B1) == value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim + 1, precision_type2(value))
            @test B1 != B2
        end
    end

    # dim(B1) == dim(B2), value(B1) != value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim, precision_type2(value) + 1)
            @test B1 != B2
        end
    end
end

@testset "!=(x, y): x, y::{AbstractBlade, Real}" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    dim = size(vectors, 1)
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- x::AbstractBlade, y::AbstractBlade

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # x::Blade, y::Scalar
            # x::Scalar, y::Blade
            B1 = Blade{precision_type1}(vectors)
            B2 = Scalar{precision_type2}(test_value)
            @test B1 != B2
            @test B2 != B1

            # x::Blade, y::Pseudoscalar
            # x::Pseudoscalar, y::Blade
            B1 = Blade{precision_type1}(vectors)
            B2 = Pseudoscalar{precision_type2}(dim, test_value)
            @test B1 != B2
            @test B2 != B1

            # x::Scalar, y::Pseudoscalar
            # x::Pseudoscalar, y::Scalar
            B1 = Scalar{precision_type1}(test_value)
            B2 = Pseudoscalar{precision_type2}(dim, test_value)
            @test B1 != B2
            @test B2 != B1
        end
    end

    # --- x::Union{Blade, Pseudoscalar}, y::Real
    #     x::Real, y::Union{Blade, Pseudoscalar}

    for precision_type in subtypes(AbstractFloat)
        # x::Blade, y::Real
        # x::Real, y::Blade
        B = Blade{precision_type}(vectors)
        @test B != test_value
        @test test_value != B

        # x::Pseudoscalar, y::Real
        # x::Real, y::Pseudoscalar
        B = Pseudoscalar{precision_type}(dim, test_value)
        @test B != test_value
        @test test_value != B
    end
end

# --- ≈(x, y)

@testset "≈(x, y): x, y::Blade" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B1) == dim(B2), grade(B1) == grade(B2), volume(B1) ≈ volume(B2)
    # basis(B1) ≈ basis(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors))
            @test B1 ≈ B2
        end
    end

    # dim(B1) != dim(B2)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors2))
            @test B1 ≉ B2
        end
    end

    # grade(B1) != grade(B2)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors2))
            @test B1 ≉ B2
        end
    end

    # volume(B1) ≉ volume(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors))
            B2 = Blade(convert(Array{precision_type2}, vectors),
                       volume=2*volume(B1))
            @test B1 ≉ B2
        end
    end

    # basis(B1) ≉ basis(B2)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 10; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            B2 = Blade(convert(Array{precision_type2}, vectors2),
                       volume=test_volume)
            @test B1 ≉ B2
        end
    end

    # B1 and B2 have opposite orientations
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            B2 = Blade(convert(Array{precision_type2}, vectors),
                       volume=-test_volume)
            @test B1 ≉ B2
        end
    end
end

@testset "≈(x, y): x, y::Scalar" begin
    # Preparations
    value = rand()
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

@testset "≈(x, y): x, y::Pseudoscalar" begin
    # Preparations
    dim = 10

    value = rand()
    value = rand() > 0.5 ? value : -value

    # dim(B1) == dim(B2), value(B1) == value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim, precision_type2(value))
            @test B1 ≈ B2
        end
    end

    # dim(B1) != dim(B2), value(B1) == value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim + 1, precision_type2(value))
            @test B1 ≉ B2
        end
    end

    # dim(B1) == dim(B2), value(B1) != value(B2)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Pseudoscalar(dim, precision_type1(value))
            B2 = Pseudoscalar(dim, precision_type2(value) + 1)
            @test B1 ≉ B2
        end
    end
end


@testset "≉(x, y): x, y::{AbstractBlade, Real}" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    dim = size(vectors, 1)
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- x::AbstractBlade, y::AbstractBlade

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # x::Blade, y::Scalar
            # x::Scalar, y::Blade
            B1 = Blade{precision_type1}(vectors)
            B2 = Scalar{precision_type2}(test_value)
            @test B1 ≉ B2
            @test B2 ≉ B1

            # x::Blade, y::Pseudoscalar
            # x::Pseudoscalar, y::Blade
            B1 = Blade{precision_type1}(vectors)
            B2 = Pseudoscalar{precision_type2}(dim, test_value)
            @test B1 ≉ B2
            @test B2 ≉ B1

            # x::Scalar, y::Pseudoscalar
            # x::Pseudoscalar, y::Scalar
            B1 = Scalar{precision_type1}(test_value)
            B2 = Pseudoscalar{precision_type2}(dim, test_value)
            @test B1 ≉ B2
            @test B2 ≉ B1
        end
    end

    # --- x::Union{Blade, Pseudoscalar}, y::Real
    #     x::Real, y::Union{Blade, Pseudoscalar}

    for precision_type in subtypes(AbstractFloat)
        # x::Blade, y::Real
        # x::Real, y::Blade
        B = Blade{precision_type}(vectors)
        @test B ≉ test_value
        @test test_value ≉ B

        # x::Pseudoscalar, y::Real
        # x::Real, y::Pseudoscalar
        B = Pseudoscalar{precision_type}(dim, test_value)
        @test B ≉ test_value
        @test test_value ≉ B
    end
end
