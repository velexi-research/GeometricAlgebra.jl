"""
Unit tests for the Pseudoscalar type.

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

@testset "Pseudoscalar: inner constructor" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # --- Pseudoscalar{T}(dim::Integer, value::Real)

    test_dim = 10
    test_value = rand() + 1  # add 1 to avoid 0
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = precision_type(test_value)

            S = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test S.dim == test_dim
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end
    end

    # --- Pseudoscalar{T}(dim::Integer, value::Integer)

    test_dim = 10
    test_value = (rand() > 0.5) ? 10 : -10

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            converted_test_value = precision_type(test_value)

            S = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test S.dim == test_dim
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = precision_type(abs(test_value))

            S = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test S.dim == test_dim
            @test S.value isa precision_type
            @test S.value == precision_type(converted_test_value)
        end

        # Bool
        S = Pseudoscalar{precision_type}(test_dim, true)
        @test S.dim == test_dim
        @test S.value isa precision_type
        @test S.value == precision_type(1)

        S = Pseudoscalar{precision_type}(test_dim, false)
        @test S == zero(Pseudoscalar{precision_type})
    end

    # --- Invalid arguments

    for precision_type in subtypes(AbstractFloat)
        # dim == 0
        @test_throws ErrorException Pseudoscalar(0, 10)

        # dim < 0
        @test_throws ErrorException Pseudoscalar(-1, 10)
    end
end

@testset "Pseudoscalar: outer constructor - basic constructors" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.
    =#

    # --- Preparations

    test_dim = 10

    # --- Pseudoscalar(dim::Integer, value::AbstractFloat)

    test_value = rand()
    test_value = (rand() > 0.5) ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        converted_test_value = precision_type(test_value)
        S = Pseudoscalar(test_dim, converted_test_value)
        @test S isa Pseudoscalar{precision_type}
    end

    # --- Pseudoscalar(dim::Integer, value::Integer)

    test_value = (rand() > 0.5) ? 10 : -10

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        S = Pseudoscalar(test_dim, value_type(test_value))
        @test S isa Pseudoscalar{Float64}
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        S = Pseudoscalar(test_dim, value_type(abs(test_value)))
        @test S isa Pseudoscalar{Float64}
    end

    # Bool
    S = Pseudoscalar(test_dim, true)
    @test S isa Pseudoscalar{Float64}

    S = Pseudoscalar(test_dim, false)
    @test S == zero(Pseudoscalar{Float64})
end

@testset "Pseudoscalar: outer constructor - copy constructor" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.

      * Test behavior of keyword arguments: `value`.
    =#

    # --- Preparations

    test_dim = 10

    test_value = rand()
    test_value = (rand() > 0.5) ? test_value : -test_value

    # --- Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)
        S = Pseudoscalar(test_dim, converted_test_value)

        # Construct a Pseudoscalar representing the same pseudoscalar as `S`
        S_copy = Pseudoscalar(S)
        @test S_copy isa Pseudoscalar{precision_type}

        # Construct a Pseudoscalar representing the same space as `S` with a
        # different value.
        S_copy = Pseudoscalar(S, value=converted_test_value + 1)
        @test S_copy isa Pseudoscalar{precision_type}
        @test value(S_copy) == converted_test_value + 1
    end
end

# --- Function tests

@testset "AbstractBlade interface: S::Pseudoscalar" begin
    # Preparations
    test_dim = 10

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        S = Pseudoscalar(test_dim, positive_test_value)
        @test dim(S) == test_dim
        @test grade(S) == test_dim
        @test basis(S) == LinearAlgebra.I
        @test volume(S) isa precision_type
        @test volume(S) == positive_test_value
        @test norm(S) isa precision_type
        @test norm(S) == positive_test_value
        @test sign(S) == 1

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Pseudoscalar(test_dim, negative_test_value)
        @test dim(S) == test_dim
        @test grade(S) == test_dim
        @test basis(S) == LinearAlgebra.I
        @test volume(S) isa precision_type
        @test volume(S) == negative_test_value
        @test norm(S) isa precision_type
        @test norm(S) == abs(negative_test_value)
        @test sign(S) == -1

        # value = Inf
        S = Pseudoscalar(test_dim, precision_type(Inf))
        @test dim(S) == test_dim
        @test grade(S) == test_dim
        @test basis(S) == LinearAlgebra.I
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(Inf)
        @test norm(S) isa precision_type
        @test norm(S) == Inf
        @test sign(S) == 1

        # value = -Inf
        S = Pseudoscalar(test_dim, precision_type(-Inf))
        @test dim(S) == test_dim
        @test grade(S) == test_dim
        @test basis(S) == LinearAlgebra.I
        @test volume(S) isa precision_type
        @test volume(S) == precision_type(-Inf)
        @test norm(S) isa precision_type
        @test norm(S) == Inf
        @test sign(S) == -1
    end
end

@testset "Pseudoscalar: Pseudoscalar interface" begin
    # Preparations
    test_dim = 10

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        S = Pseudoscalar(test_dim, positive_test_value)
        @test value(S) isa precision_type
        @test value(S) == positive_test_value

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        S = Pseudoscalar(test_dim, negative_test_value)
        @test value(S) isa precision_type
        @test value(S) == negative_test_value

        # value = 0
        S = Pseudoscalar(test_dim, precision_type(0))
        @test value(S) isa precision_type
        @test value(S) == 0

        # value = Inf
        S = Pseudoscalar(test_dim, precision_type(Inf))
        @test value(S) isa precision_type
        @test value(S) == precision_type(Inf)

        # value = -Inf
        S = Pseudoscalar(test_dim, precision_type(-Inf))
        @test value(S) isa precision_type
        @test value(S) == precision_type(-Inf)
    end
end

@testset "convert(S): S::Pseudoscalar" begin
    # Preparations
    test_dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # Tests
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            converted_test_value = precision_type_src(test_value)
            S = Pseudoscalar{precision_type_src}(test_dim, converted_test_value)

            # Exercise functionality and check results
            S_converted = convert(Pseudoscalar{precision_type_converted}, S)

            @test S_converted isa Pseudoscalar{precision_type_converted}

            if precision_type_src == precision_type_converted
                @test S_converted === S
            else
                @test S_converted !== S
                @test S_converted ≈ S
            end
        end
    end
end
