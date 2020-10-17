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
using Test
import Random.rand

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Constructor tests

@testset "Scalar: inner constructor" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # --- Scalar{T}(value::Real)

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)

            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end
    end

    # --- Scalar{T}(value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            converted_test_value = value_type(test_value)

            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = value_type(abs(test_value))

            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end

        # Bool
        S = Scalar{precision_type}(true)
        @test S.value isa precision_type
        @test S.value == precision_type(1)

        S = Scalar{precision_type}(false)
        @test S.value isa precision_type
        @test S.value == precision_type(0)
    end
end

@testset "Scalar: outer constructor - basic constructors" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.
    =#

    # --- Scalar(value::AbstractFloat)

    test_value = rand()
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        converted_test_value = precision_type(test_value)
        S = Scalar(converted_test_value)
        @test S isa Scalar{precision_type}
    end

    # --- Scalar(value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        S = Scalar(value_type(test_value))
        @test S isa Scalar{Float64}
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        S = Scalar(value_type(abs(test_value)))
        @test S isa Scalar{Float64}
    end

    # Bool
    S = Scalar(true)
    @test S isa Scalar{Float64}

    S = Scalar(false)
    @test S isa Scalar{Float64}
end

@testset "Scalar: outer constructor - copy constructor" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.

      * Test behavior of keyword arguments: `value`.
    =#

    # --- Preparations

    test_value = rand()
    test_value = (rand() > 0.5) ? test_value : -test_value

    # --- Scalar(S::Scalar, value::Real)

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)
        S = Scalar(converted_test_value)

        # Construct a Scalar representing the same scalar as `S`
        S_copy = Scalar(S)
        @test S_copy isa Scalar{precision_type}

        # Construct a Scalar representing the same space as `S` with a
        # different value.
        S_copy = Scalar(S, value=converted_test_value + 1)
        @test S_copy isa Scalar{precision_type}
        @test value(S_copy) == converted_test_value + 1
    end
end

# --- Function tests

@testset "AbstractBlade interface: S::Scalar" begin
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
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == positive_test_value
        @test norm(S) isa precision_type
        @test norm(S) == positive_test_value
        @test sign(S) == 1

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Scalar(negative_test_value)
        @test dim(S) == 0
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == negative_test_value
        @test norm(S) isa precision_type
        @test norm(S) == abs(negative_test_value)
        @test sign(S) == -1

        # value = 0
        S = Scalar(precision_type(0))
        @test dim(S) == 0
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == 0
        @test norm(S) isa precision_type
        @test norm(S) == 0
        @test sign(S) == 0

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test dim(S) == 0
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(Inf)
        @test norm(S) isa precision_type
        @test norm(S) == Inf
        @test sign(S) == 1

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test dim(S) == 0
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(-Inf)
        @test norm(S) isa precision_type
        @test norm(S) == Inf
        @test sign(S) == -1
    end
end

@testset "Scalar: Scalar interface" begin
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
        @test value(S) isa precision_type
        @test value(S) == positive_test_value

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Scalar(negative_test_value)
        @test value(S) isa precision_type
        @test value(S) == negative_test_value

        # value = 0
        S = Scalar(precision_type(0))
        @test value(S) isa precision_type
        @test value(S) == 0

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test value(S) isa precision_type
        @test value(S) == precision_type(Inf)

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test value(S) isa precision_type
        @test value(S) == precision_type(-Inf)
    end
end

@testset "convert(S): S::Scalar" begin
    # Preparations
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # Tests
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            converted_test_value = precision_type_src(test_value)
            S = Scalar{precision_type_src}(converted_test_value)

            # Exercise functionality and check results
            S_converted = convert(Scalar{precision_type_converted}, S)

            @test S_converted isa Scalar{precision_type_converted}

            if precision_type_src == precision_type_converted
                @test S_converted === S
            else
                @test S_converted !== S
                @test S_converted ≈ S
            end
        end
    end
end
