"""
Unit tests for the Scalar type.

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
import Random.rand

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Constructor tests

@testset "Scalar: inner constructor tests" begin
    # Notes
    # -----
    # * Test value of constructed instance

    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = (rand() > 0.5) ? value : -value

    # Scalar{T}(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_value = precision_type(value)

        # Default atol
        S = Scalar{precision_type}(converted_value)
        @test S.value == converted_value
        @test typeof(S.value) === precision_type

        # abs(value) > atol
        S = Scalar{precision_type}(converted_value,
                                   atol=abs(converted_value)-1)
        @test S.value == converted_value
        @test typeof(S.value) === precision_type

        # abs(value) == atol
        S = Scalar{precision_type}(converted_value, atol=abs(converted_value))
        @test S.value == converted_value
        @test typeof(S.value) === precision_type

        # abs(value) < atol
        S = Scalar{precision_type}(converted_value,
                                   atol=abs(converted_value)+1)
        @test S == Zero(precision_type)
    end
end

@testset "Scalar: outer constructor tests" begin
    # Notes
    # -----
    # * Test type of constructed instances. Correct construction of instances
    #   is tested by the inner constructor tests.
    #
    # * Test behavior of `atol` argument.

    # Preparations
    value = -3

    # --- Scalar(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_value = precision_type(value)

        # Default atol
        S = Scalar(converted_value)
        @test typeof(S) === Scalar{precision_type}

        # abs(value) > atol
        S = Scalar(converted_value, atol=abs(value)-1)
        @test typeof(S) === Scalar{precision_type}

        # abs(value) == atol
        S = Scalar(converted_value, atol=abs(value))
        @test typeof(S) === Scalar{precision_type}

        # abs(value) < atol
        S = Scalar(converted_value, atol=abs(value)+1)
        @test S === Zero(precision_type)
    end

    # --- Scalar{T}(value::S;
    #               atol::Real=eps(T)) where {T<:AbstractFloat,
    #                                         S<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        for value_precision_type in subtypes(AbstractFloat)
            # Preparations
            converted_value = value_precision_type(value)

            # Default atol
            S = Scalar{precision_type}(converted_value)
            @test typeof(S) === Scalar{precision_type}

            # abs(value) > atol
            S = Scalar{precision_type}(converted_value, atol=abs(value)-1)
            @test typeof(S) === Scalar{precision_type}

            # abs(value) == atol
            S = Scalar{precision_type}(converted_value, atol=abs(value))
            @test typeof(S) === Scalar{precision_type}

            # abs(value) < atol
            S = Scalar{precision_type}(converted_value, atol=abs(value)+1)
            @test S === Zero(precision_type)
        end
    end

    # --- Scalar(value::Integer)

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        # value is nonzero
        S = Scalar(convert(value_type, value))
        @test typeof(S) === Scalar{Float64}

        # value is zero
        S = Scalar(convert(value_type, 0))
        @test S === Zero(Float64)
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        # value is nonzero
        S = Scalar(convert(value_type, abs(value)))
        @test typeof(S) === Scalar{Float64}

        # value is zero
        S = Scalar(convert(value_type, 0))
        @test S === Zero(Float64)
    end

    # Bool
    S = Scalar(true)
    @test typeof(S) === Scalar{Float64}

    S = Scalar(false)
    @test S === Zero(Float64)

    # --- Scalar{T}(value::Integer) where {T<:AbstractFloat}

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        for precision_type in subtypes(AbstractFloat)
            # value is nonzero
            S = Scalar{precision_type}(convert(value_type, value))
            @test typeof(S) === Scalar{precision_type}

            # value is zero
            S = Scalar{precision_type}(convert(value_type, 0))
            @test S === Zero(precision_type)
        end
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        for precision_type in subtypes(AbstractFloat)
            # value is nonzero
            S = Scalar{precision_type}(convert(value_type, abs(value)))
            @test typeof(S) === Scalar{precision_type}

            # value is zero
            S = Scalar{precision_type}(convert(value_type, 0))
            @test S === Zero(precision_type)
        end
    end

    # Bool
    for precision_type in subtypes(AbstractFloat)
        # value is nonzero
        S = Scalar{precision_type}(true)
        @test typeof(S) === Scalar{precision_type}

        # value is zero
        value = 0
        S = Scalar{precision_type}(false)
        @test S === Zero(precision_type)
    end
end


# --- Function tests

@testset "Scalar: function tests" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_value = precision_type(value)

        # value > 0
        B = Scalar(converted_value)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == abs(converted_value)
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(1 / converted_value)

        # value < 0
        negative_value = -(abs(converted_value))
        B = Scalar(negative_value)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == abs(negative_value)
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(1 / negative_value)

        # value = 0
        B = Scalar(precision_type(0))
        @test B === Zero(precision_type)
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == 0
        @test basis(B) === nothing
        @test inverse(B) == Scalar{precision_type}(Inf)

        # value = Inf
        B = Scalar(precision_type(Inf))
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == Inf
        @test basis(B) === nothing
        @test inverse(B) === Zero(precision_type)

        # value = -Inf
        B = Scalar(precision_type(-Inf))
        @test dim(B) == 0
        @test grade(B) == 0
        @test norm(B) == Inf
        @test basis(B) === nothing
        @test inverse(B) === Zero(precision_type)
    end
end

@testset "Scalar: comparison operation tests" begin
    # Preparations
    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value

    float64_or_bigfloat = (Float64, BigFloat)

    # :(==)
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        B = Scalar(converted_value)
        @test B == converted_value
        @test converted_value == B
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
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

    # :(≈)
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(value)
        B = Scalar(converted_value)
        @test B ≈ converted_value
        @test converted_value ≈ B
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B1 = Scalar(precision_type1(value))
            B2 = Scalar(precision_type2(value))
            @test B1 ≈ B2
        end
    end
end
