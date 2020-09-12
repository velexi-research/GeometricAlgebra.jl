"""
Unit tests for the Scalar type.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test
import LinearAlgebra
import Random.rand

# GeometricAlgebra.jl
using GeometricAlgebra
import GeometricAlgebra.≈


# --- Constructor tests

@testset "Scalar: constructor tests" begin
    # Scalar(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        # Default atol
        value = 3
        S = Scalar(convert(precision_type, value))

        @test S.value == value
        @test typeof(S.value) === precision_type
        @test S == value
        @test value == S

        # atol < abs(value)
        value = -3
        S = Scalar(convert(precision_type, value), atol=abs(value)-1)

        @test S.value == value
        @test typeof(S.value) === precision_type
        @test S == value
        @test value == S

        # atol > abs(value)
        value = -3
        S = Scalar(convert(precision_type, value), atol=abs(value)+1)
        @test S == Zero(precision_type)
    end

    # Scalar{T}(value::S;
    #           atol::Real=eps(T)) where {T<:AbstractFloat, S<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        for value_precision_type in subtypes(AbstractFloat)
            # Default atol
            value = 3
            S = Scalar{precision_type}(convert(value_precision_type, value))

            @test S.value == value
            @test typeof(S.value) === precision_type
            @test S == value
            @test value == S

            # atol < abs(value)
            value = -3
            S = Scalar{precision_type}(convert(value_precision_type, value),
                                       atol=abs(value)-1)

            @test S.value == value
            @test typeof(S.value) === precision_type
            @test S == value
            @test value == S

            # atol > abs(value)
            value = -3
            S = Scalar{precision_type}(convert(value_precision_type, value),
                                       atol=abs(value)+1)
            @test S == Zero(precision_type)
        end
    end

    # Scalar(value::Integer)
    for value_type in vcat(subtypes(Unsigned), subtypes(Signed), Bool)
        # value is nonzero
        value = (value_type != Bool) ? 3 : 1
        S = Scalar(convert(value_type, value))

        @test S.value == value
        @test typeof(S.value) === Float64
        @test S == value
        @test value == S

        # value is zero
        value = 0
        S = Scalar(convert(value_type, value))

        @test S === Zero(Float64)
    end

    # Scalar{T}(value::Integer) where {T<:AbstractFloat}
    for value_type in vcat(subtypes(Unsigned), subtypes(Signed), Bool)
        for precision_type in subtypes(AbstractFloat)
            # value is nonzero
            value = (value_type != Bool) ? 3 : 1
            S = Scalar{precision_type}(convert(value_type, value))

            @test S.value == value
            @test typeof(S.value) === precision_type
            @test S == value
            @test value == S

            # value is zero
            value = 0
            S = Scalar{precision_type}(convert(value_type, value))

            @test S === Zero(precision_type)
        end
    end
end


# --- Function tests

@testset "Scalar: function tests" begin
    for precision_type in subtypes(AbstractFloat)
        # value > 0
        value = rand(precision_type) + 1  # avoid 0
        B = Scalar(value)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == value
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(1 / value)

        # value < 0
        value = -(rand(precision_type) + 1)  # avoid 0
        B = Scalar(precision_type(value))
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == abs(value)
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(1 / precision_type(value))

        # value = 0
        value = precision_type(0)
        B = Scalar(value)
        @test B === Zero(precision_type)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == 0
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(Inf)

        # value = Inf
        value = precision_type(Inf)
        B = Scalar(value)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == Inf
        @test basis(B) === nothing
        @test inverse(B) === Zero(precision_type)

        # value = -Inf
        value = -precision_type(Inf)
        B = Scalar(value)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == Inf
        @test basis(B) === nothing
        @test inverse(B) === Zero(precision_type)
    end
end

@testset "Scalar: comparison operation tests" begin
    float64_or_bigfloat = (Float64, BigFloat)

    # :(==)
    for precision_type in subtypes(AbstractFloat)
        value = rand(precision_type) + 1  # avoid 0
        B = Scalar(value)
        @test B == value
        @test value == B
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            value = rand() + 1  # avoid 0
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
        end
    end

    # isapprox()
    for precision_type in subtypes(AbstractFloat)
        value = rand(precision_type) + 1  # avoid 0
        B = Scalar(value)
        @test B ≈ value
        @test value ≈ B
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            value = rand() + 1  # avoid 0
            B1 = Scalar(precision_type1(value))
            B2 = Scalar(precision_type2(value))
            @test B1 ≈ B2
        end
    end
end
