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

# --- File inclusions

# Test utilities
include("test_utils.jl")

# --- Constructor tests

@testset "Scalar: inner constructor" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # --- Scalar{T}(value::AbstractFloat)

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)

            # value != 0, value != 1
            S = Scalar{precision_type}(converted_test_value)
            @test S isa Scalar{precision_type}
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            @test S isa Zero{precision_type}

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            @test S isa One{precision_type}
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
            S isa Scalar{precision_type}
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            S isa Zero{precision_type}

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            S isa One{precision_type}
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = value_type(abs(test_value))

            # value != 0, value != 1
            S = Scalar{precision_type}(converted_test_value)
            S isa Scalar{precision_type}
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)

            # value = 0
            S = Scalar{precision_type}(value_type(0))
            @test S isa Zero{precision_type}

            # value = 1
            S = Scalar{precision_type}(value_type(1))
            S isa One{precision_type}
        end

        # Bool
        S = Scalar{precision_type}(true)
        @test S isa One{precision_type}

        S = Scalar{precision_type}(false)
        @test S isa Zero{precision_type}
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

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

    for precision_type in subtypes(AbstractFloat)
        # value != 0, value != 1
        converted_test_value = precision_type(test_value)
        S = Scalar(converted_test_value)
        @test S isa Scalar{precision_type}
        @test S.value == precision_type(converted_test_value)

        # value = 0
        converted_test_value = precision_type(0)
        S = Scalar(converted_test_value)
        @test S isa Zero{precision_type}

        # value = 1
        converted_test_value = precision_type(1)
        S = Scalar(converted_test_value)
        @test S isa One{precision_type}
    end

    # --- Scalar(value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        # value != 0, value != 1
        S = Scalar(value_type(test_value))
        @test S isa Scalar{Float64}
        @test S.value == test_value

        # value = 0
        converted_test_value = value_type(0)
        S = Scalar(converted_test_value)
        @test S isa Zero{Float64}

        # value = 1
        converted_test_value = value_type(1)
        S = Scalar(converted_test_value)
        @test S isa One{Float64}
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        # value != 0, value != 1
        S = Scalar(value_type(abs(test_value)))
        @test S isa Scalar{Float64}
        @test S.value == abs(test_value)

        # value = 0
        converted_test_value = value_type(0)
        S = Scalar(converted_test_value)
        @test S isa Zero{Float64}

        # value = 1
        converted_test_value = value_type(1)
        S = Scalar(converted_test_value)
        @test S isa One{Float64}
    end

    # Bool
    S = Scalar(true)
    @test S isa One{Float64}

    S = Scalar(false)
    @test S isa Zero{Float64}
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
        @test iszero(S_copy)

        # value = 1
        S = Scalar(precision_type(1))
        S_copy = Scalar(S)
        @test isone(S_copy)
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

# --- Test attribute methods

@testset "Scalar: AbstractMultivector attribute functions" begin
    # --- Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

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
        @test norm(S) isa precision_type
        @test norm(S) == abs(converted_test_value)

        @test S[0] == [S]
        @test S[1] == []

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

@testset "Scalar: AbstractBlade attribute functions" begin
    # --- Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

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

@testset "Scalar: AbstractScalar attribute functions" begin
    # Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

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

# --- Test comparison operations

@testset "Scalar: ==(B, C)" begin
    # Preparations

    test_value = get_random_value()
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

@testset "Scalar: ≈(B, C)" begin
    # Preparations

    test_value = get_random_value()

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

# --- Tests for AbstractMultivector interface functions

@testset "Scalar: -(B)" begin
    # Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        minus_B = -B
        @test minus_B isa Scalar{precision_type}
        @test minus_B == Scalar{precision_type}(-test_value)
    end
end

@testset "Scalar: reverse(B)" begin
    # --- Preparations

    test_value = get_random_value(1) # add 1 to avoid 0

    # --- Exercise functionality and check results

    # value != 0
    B = Scalar(test_value)
    reverse_B = reverse(B)
    @test reverse_B === B
    @test B * reverse_B ≈ norm(B)^2

    # value = Inf
    B = Scalar(Inf)
    @test reverse(B) === B

    # value = -Inf
    B = Scalar(-Inf)
    @test reverse(B) === B
end

@testset "Scalar: dual(B)" begin
    # Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)

            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(test_value)) :
                Pseudoscalar{precision_type}(test_dim,
                                             precision_type(-test_value))
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_result
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "Scalar: convert(S)" begin
    # Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

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

# --- Tests for AbstractBlade interface functions

@testset "Scalar: reciprocal(B)" begin
    # Preparations

    test_value = get_random_value(2) # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)

        # value > 0
        B = Scalar(abs(converted_value))
        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar(1 / abs(converted_value))
        @test B * reciprocal_B ≈ 1

        # value < 0
        negative_value = -(abs(converted_value))
        B = Scalar(negative_value)
        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar(1 / negative_value)
        @test B * reciprocal_B ≈ 1

        # value = Inf
        B = Scalar(precision_type(Inf))
        @test iszero(reciprocal(B))

        # value = -Inf
        B = Scalar(precision_type(-Inf))
        @test iszero(reciprocal(B))
    end
end
