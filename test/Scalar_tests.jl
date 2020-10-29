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

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)

            # value != 0, value != 1
            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            @test S === Zero{precision_type}()

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            @test S === One{precision_type}()
        end
    end

    # --- Scalar{T}(value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            converted_test_value = value_type(test_value)

            # value != 0, value != 1
            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            @test S === Zero{precision_type}()

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            @test S === One{precision_type}()
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = value_type(abs(test_value))

            # value != 0, value != 1
            S = Scalar{precision_type}(converted_test_value)
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            @test S === Zero{precision_type}()

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            @test S === One{precision_type}()
        end

        # Bool
        S = Scalar{precision_type}(true)
        @test S === One{precision_type}()

        S = Scalar{precision_type}(false)
        @test S === Zero{precision_type}()
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

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        # value != 0, value != 1
        converted_test_value = precision_type(test_value)
        S = Scalar(converted_test_value)
        @test S isa Scalar{precision_type}

        # value = 0
        converted_test_value = precision_type(0)
        S = Scalar(converted_test_value)
        @test S === Zero{precision_type}()

        # value = 1
        converted_test_value = precision_type(1)
        S = Scalar(converted_test_value)
        @test S === One{precision_type}()
    end

    # --- Scalar(value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        # value != 0, value != 1
        S = Scalar(value_type(test_value))
        @test S isa Scalar{Float64}

        # value = 0
        converted_test_value = value_type(0)
        S = Scalar(converted_test_value)
        @test S === Zero{Float64}()

        # value = 1
        converted_test_value = value_type(1)
        S = Scalar(converted_test_value)
        @test S === One{Float64}()
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        # value != 0, value != 1
        S = Scalar(value_type(abs(test_value)))
        @test S isa Scalar{Float64}

        # value = 0
        converted_test_value = value_type(0)
        S = Scalar(converted_test_value)
        @test S === Zero{Float64}()

        # value = 1
        converted_test_value = value_type(1)
        S = Scalar(converted_test_value)
        @test S === One{Float64}()
    end

    # Bool
    S = Scalar(true)
    @test S === One{Float64}()

    S = Scalar(false)
    @test S === Zero{Float64}()
end

#= DEPRECATED
@testset "Scalar: outer constructor - copy constructor" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.

      * Test behavior of keyword arguments: `value`.
    =#

    # --- Preparations

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # --- Scalar(S::Scalar)

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value != 0, value != 1
        S = Scalar(converted_test_value)
        S_copy = Scalar(S)
        @test S_copy isa Scalar{precision_type}

        # value = 0
        S = Scalar(precision_type(0))
        S_copy = Scalar(S)
        @test S_copy === Zero{precision_type}()

        # value = 1
        S = Scalar(precision_type(1))
        S_copy = Scalar(S)
        @test S_copy === One{precision_type}()
    end

    # --- Scalar(S::Scalar; value::Real)

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value != 0, value != 1
        S = Scalar(converted_test_value)
        S_copy = Scalar(S, value=converted_test_value + 1)
        @test S_copy isa Scalar{precision_type}
        @test value(S_copy) == converted_test_value + 1

        # value = 0
        S = Scalar(precision_type(0))
        S_copy = Scalar(S, value=converted_test_value + 1)
        @test S_copy isa Scalar{precision_type}
        @test value(S_copy) == converted_test_value + 1

        # value = 1
        S = Scalar(precision_type(1))
        S_copy = Scalar(S, value=converted_test_value + 1)
        @test S_copy isa Scalar{precision_type}
        @test value(S_copy) == converted_test_value + 1
    end
end
=#

# --- Function tests

@testset "Scalar: AbstractMultivector interface functions" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # --- Tests

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        S = Scalar{precision_type}(positive_test_value)
        @test dim(S) == 0
        @test grades(S) == [0]
        @test blades(S) == [S]
        @test S[0] == [S]
        @test S[1] == []
        @test norm(S) isa precision_type
        @test norm(S) == abs(converted_test_value)

        # value < 0
        negative_test_value = converted_test_value < 0 ?
            converted_test_value : -converted_test_value
        S = Scalar{precision_type}(negative_test_value)
        @test norm(S) isa precision_type
        @test norm(S) == abs(converted_test_value)

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test norm(S) isa precision_type
        @test norm(S) == Inf

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test norm(S) isa precision_type
        @test norm(S) == Inf
    end
end

@testset "Scalar: AbstractBlade interface functions" begin
    # --- Preparations

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # --- Tests

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        S = Scalar(positive_test_value)
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == positive_test_value
        @test sign(S) == 1

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Scalar(negative_test_value)
        @test volume(S) isa precision_type
        @test volume(S) == negative_test_value
        @test sign(S) == -1

        # value = Inf
        S = Scalar(precision_type(Inf))
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(Inf)
        @test sign(S) == 1

        # value = -Inf
        S = Scalar(precision_type(-Inf))
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(-Inf)
        @test sign(S) == -1
    end
end

@testset "Scalar: AbstractScalar interface functions" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Exercise functionality and check results
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

@testset "Scalar: convert(S)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Exercise functionality and check results
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            converted_test_value = precision_type_src(test_value)
            S = Scalar{precision_type_src}(converted_test_value)
            S_converted = convert(AbstractScalar{precision_type_converted}, S)
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
