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
    test_value = rand() + 1  # add 1 to avoid 0
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Scalar{T}(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # abs(value) > default atol
        S = Scalar{precision_type}(converted_test_value)
        @test S.value == converted_test_value
        @test S.value isa precision_type

        # abs(value) < default atol
        S = Scalar{precision_type}(blade_atol(precision_type) / 2)
        @test S === Zero{precision_type}()

        # abs(value) > atol
        S = Scalar{precision_type}(converted_test_value,
                                   atol=abs(converted_test_value) - 1)
        @test S.value == converted_test_value
        @test S.value isa precision_type

        # abs(value) == atol
        S = Scalar{precision_type}(converted_test_value,
                                   atol=abs(converted_test_value))
        @test S.value == converted_test_value
        @test S.value isa precision_type

        # abs(value) < atol
        S = Scalar{precision_type}(converted_test_value,
                                   atol=abs(converted_test_value) + 1)
        @test S === Zero{precision_type}()
    end
end

@testset "Scalar: outer constructor tests" begin
    # Notes
    # -----
    # * Test type of constructed instances. Correct construction of instances
    #   is tested by the inner constructor tests.
    #
    # * Test behavior of `atol` argument.

    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = (rand() > 0.5) ? test_value : -test_value

    int_test_value = (rand() > 0.5) ? 10 : -10

    # --- Scalar(value::T; atol::Real=eps(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # abs(value) > default atol
        S = Scalar(converted_test_value)
        @test S isa Scalar{precision_type}

        # abs(value) < default atol
        S = Scalar(precision_type(blade_atol(precision_type) / 2))
        @test S === Zero{precision_type}()

        # abs(value) > atol
        S = Scalar(converted_test_value, atol=abs(converted_test_value) - 1)
        @test S isa Scalar{precision_type}

        # abs(value) == atol
        S = Scalar(converted_test_value, atol=abs(converted_test_value))
        @test S isa Scalar{precision_type}

        # abs(value) < atol
        S = Scalar(converted_test_value, atol=abs(converted_test_value) + 1)
        @test S === Zero{precision_type}()
    end

    # --- Scalar{T}(value::AbstractFloat;
    #               atol::Real=eps(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            # Preparations
            converted_test_value = value_type(test_value)

            # abs(value) > default atol
            S = Scalar{precision_type}(converted_test_value)
            @test S isa Scalar{precision_type}

            # abs(value) < default atol
            S = Scalar{precision_type}(
                value_type(blade_atol(precision_type) / 2))
            @test S === Zero{precision_type}()

            # abs(value) > atol
            S = Scalar{precision_type}(
                converted_test_value,
                atol=abs(precision_type(converted_test_value)) - 1)
            @test S isa Scalar{precision_type}

            # abs(value) == atol
            S = Scalar{precision_type}(
                converted_test_value,
                atol=abs(precision_type(converted_test_value)))
            @test S isa Scalar{precision_type}

            # abs(value) < atol
            S = Scalar{precision_type}(
                converted_test_value,
                atol=abs(precision_type(converted_test_value)) + 1)
            @test S === Zero{precision_type}()
        end
    end

    # --- Scalar(value::Integer)

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        # value is nonzero
        S = Scalar(convert(value_type, int_test_value))
        @test S isa Scalar{Float64}

        # value is zero
        S = Scalar(convert(value_type, 0))
        @test S === Zero{Float64}()
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        # value is nonzero
        S = Scalar(convert(value_type, abs(int_test_value)))
        @test S isa Scalar{Float64}

        # value is zero
        S = Scalar(convert(value_type, 0))
        @test S === Zero{Float64}()
    end

    # Bool
    S = Scalar(true)
    @test S isa Scalar{Float64}

    S = Scalar(false)
    @test S === Zero{Float64}()

    # --- Scalar{T}(value::Integer) where {T<:AbstractFloat}

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        for precision_type in subtypes(AbstractFloat)
            # value is nonzero
            S = Scalar{precision_type}(convert(value_type, int_test_value))
            @test S isa Scalar{precision_type}

            # value is zero
            S = Scalar{precision_type}(convert(value_type, 0))
            @test S === Zero{precision_type}()
        end
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        for precision_type in subtypes(AbstractFloat)
            # value is nonzero
            S = Scalar{precision_type}(convert(value_type, abs(int_test_value)))
            @test S isa Scalar{precision_type}

            # value is zero
            S = Scalar{precision_type}(convert(value_type, 0))
            @test S === Zero{precision_type}()
        end
    end

    # Bool
    for precision_type in subtypes(AbstractFloat)
        # value is nonzero
        S = Scalar{precision_type}(true)
        @test S isa Scalar{precision_type}

        # value is zero
        value = 0
        S = Scalar{precision_type}(false)
        @test S === Zero{precision_type}()
    end
end


# --- Function tests

@testset "Scalar: basic function tests" begin
    # Preparations
    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        S = Scalar(positive_test_value)
        @test dim(S) == 0
        @test grade(S) == 0
        @test norm(S) == positive_test_value
        @test basis(S) === nothing
        @test value(S) isa precision_type
        @test value(S) == positive_test_value

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Scalar(negative_test_value)
        @test dim(S) == 0
        @test grade(S) == 0
        @test norm(S) == abs(negative_test_value)
        @test basis(S) === nothing
        @test value(S) isa precision_type
        @test value(S) == negative_test_value

        # value = 0
        S = Scalar(precision_type(0))
        @test S === Zero{precision_type}()
        @test dim(S) == 0
        @test grade(S) == 0
        @test norm(S) == 0
        @test basis(S) === nothing
        @test value(S) isa precision_type
        @test value(S) == 0

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test dim(S) == 0
        @test grade(S) == 0
        @test norm(S) == Inf
        @test basis(S) === nothing
        @test value(S) isa precision_type
        @test value(S) == precision_type(Inf)

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test dim(S) == 0
        @test grade(S) == 0
        @test norm(S) == Inf
        @test basis(S) === nothing
        @test value(S) isa precision_type
        @test value(S) == precision_type(-Inf)
    end
end
