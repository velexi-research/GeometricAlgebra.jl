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

# GeometricAlgebra.jl
using GeometricAlgebra


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
    # Nonzero value
    value = 10
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == value
    @test basis(B) === nothing
    @test inverse(B) == Scalar(1 / value)

    # Inf
    value = Inf
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == Inf
    @test basis(B) === nothing
    @test inverse(B) === Zero()

    # -Inf
    value = -Inf
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == Inf
    @test basis(B) === nothing
    @test inverse(B) === Zero()
end

@testset "Scalar: comparison operation tests" begin
end
