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


# --- ==(B, C)

@testset "==(B, C): B, C::Blade" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) == volume(C)
    # basis(B) == basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            if precision_type1 == precision_type2
                @test B == C
            else
                @test B != C
            end
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # volume(B) != volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                       volume=2*volume(B))
            @test B != C
        end
    end

    # basis(B) != basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 5; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                       volume=test_volume)
            @test B != C
        end
    end
end

@testset "==(B, C): B, C::Scalar" begin
    # Preparations
    test_dim = 3

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    float64_or_bigfloat = (Float64, BigFloat)

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(B, C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            if precision_type1 == precision_type2
                @test B == C
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B == C
            else
                @test B != C
            end

            # value(B) != value(C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(2 * test_value))
            @test B != C
        end
    end

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    for precision_type in subtypes(AbstractFloat)
        # B::Scalar, C::AbstractFloat
        # B::AbstractFloat, C::Scalar
        for value_type in subtypes(AbstractFloat)
            B = Scalar(precision_type(test_value))

            # ==(B, C)
            if precision_type == value_type
                @test B == value_type(test_value)
                @test value_type(test_value) == B
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test B == value_type(test_value)
                @test value_type(test_value) == B
            else
                @test B != value_type(test_value)
                @test value_type(test_value) != B
            end

            # !=(x,y)
            @test B != value_type(2 * test_value)
            @test value_type(2 * test_value) != B
        end

        # B::Scalar, C::Integer
        # B::Integer, C::Scalar
        int_value::Int = 3
        for value_type in subtypes(Signed)
            B = Scalar(precision_type(int_value))

            # ==(B, C)
            @test B == value_type(int_value)
            @test value_type(int_value) == B

            # !=(x,y)
            @test B != value_type(2 * int_value)
            @test value_type(2 * int_value) != B
        end

        for value_type in subtypes(Unsigned)
            B = Scalar(precision_type(int_value))

            # ==(B, C)
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

@testset "==(B, C): B, C::Pseudoscalar" begin
    # Preparations
    dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    float64_or_bigfloat = (Float64, BigFloat)

    # dim(B) == dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value))
            if precision_type1 == precision_type2
                @test B == C
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B == C
            else
                @test B != C
            end
        end
    end

    # dim(B) != dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim + 1, precision_type2(test_value))
            @test B != C
        end
    end

    # dim(B) == dim(C), value(B) != value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value) + 1)
            @test B != C
        end
    end
end

@testset "!=(B, C): B, C::{AbstractBlade, Real}" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- B::AbstractBlade, C::AbstractBlade

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # B::Blade, C::Scalar
            # B::Scalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Scalar{precision_type2}(test_value)
            @test B != C
            @test C != B

            # B::Blade, C::Pseudoscalar
            # B::Pseudoscalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
            @test C != B

            # B::Scalar, C::Pseudoscalar
            # B::Pseudoscalar, C::Scalar
            B = Scalar{precision_type1}(test_value)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
            @test C != B
        end
    end

    # --- B::Union{Blade, Pseudoscalar}, C::Real
    #     B::Real, C::Union{Blade, Pseudoscalar}

    for precision_type in subtypes(AbstractFloat)
        # B::Blade, C::Real
        # B::Real, C::Blade
        B = Blade{precision_type}(vectors)
        @test B != test_value
        @test test_value != B

        # B::Pseudoscalar, C::Real
        # B::Real, C::Pseudoscalar
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != test_value
        @test test_value != B
    end
end

# --- ≈(B, C)

@testset "≈(B, C): B, C::Blade" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) ≈ volume(C)
    # basis(B) ≈ basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            @test B ≈ C
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # volume(B) ≉ volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                       volume=2*volume(B))
            @test B ≉ C
        end
    end

    # basis(B) ≉ basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 10; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                       volume=test_volume)
            @test B ≉ C
        end
    end

    # B and C have opposite orientations
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                       volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors),
                       volume=-test_volume)
            @test B ≉ C
        end
    end
end

@testset "≈(B, C): B, C::Scalar" begin
    # Preparations
    test_dim = 5

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            @test B ≈ C
        end
    end

    # B::Scalar, C::Real
    # B::Real, C::Scalar
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)
        B = Scalar(converted_value)
        @test B ≈ converted_value
        @test converted_value ≈ B
    end
end

@testset "≈(B, C): B, C::Pseudoscalar" begin
    # Preparations
    dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # dim(B) == dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value))
            @test B ≈ C
        end
    end

    # dim(B) != dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim + 1, precision_type2(test_value))
            @test B ≉ C
        end
    end

    # dim(B) == dim(C), value(B) != value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value) + 1)
            @test B ≉ C
        end
    end
end


@testset "≉(B, C): B, C::{AbstractBlade, Real}" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- B::AbstractBlade, C::AbstractBlade

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # B::Blade, C::Scalar
            # B::Scalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Scalar{precision_type2}(test_value)
            @test B ≉ C
            @test C ≉ B

            # B::Blade, C::Pseudoscalar
            # B::Pseudoscalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B ≉ C
            @test C ≉ B

            # B::Scalar, C::Pseudoscalar
            # B::Pseudoscalar, C::Scalar
            B = Scalar{precision_type1}(test_value)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B ≉ C
            @test C ≉ B
        end
    end

    # --- B::Union{Blade, Pseudoscalar}, C::Real
    #     B::Real, C::Union{Blade, Pseudoscalar}

    for precision_type in subtypes(AbstractFloat)
        # B::Blade, C::Real
        # B::Real, C::Blade
        B = Blade{precision_type}(vectors)
        @test B ≉ test_value
        @test test_value ≉ B

        # B::Pseudoscalar, C::Real
        # B::Real, C::Pseudoscalar
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ test_value
        @test test_value ≉ B
    end
end
