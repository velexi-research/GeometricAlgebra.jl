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

@testset "Pseudoscalar: inner constructor tests" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # Preparations
    test_dim = 10

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Pseudoscalar{T}(dim::Integer, value::T) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        S = Pseudoscalar{precision_type}(test_dim, converted_test_value)
        @test S.dim == test_dim
        @test S.value isa precision_type
        @test S.value == converted_test_value
    end
end

@testset "Pseudoscalar: outer constructor tests" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.
    =#

    # --- Preparations

    test_dim = 10

    test_value = rand()
    test_value = (rand() > 0.5) ? test_value : -test_value

    int_test_value = (rand() > 0.5) ? 10 : -10

    # --- Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        converted_test_value = precision_type(test_value)
        S = Pseudoscalar(test_dim, converted_test_value)
        @test S isa Pseudoscalar{precision_type}
    end

    # --- Pseudoscalar{T}(dim::Integer, value::AbstractFloat)
    #         where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)
            # Note: precision_type == value_type is covered by inner constructor

            if precision_type != value_type
                S = Pseudoscalar{precision_type}(test_dim, converted_test_value)
                @test S isa Pseudoscalar{precision_type}
            end
        end
    end

    # --- Pseudoscalar(dim::Integer, value::Integer)

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        S = Pseudoscalar(test_dim, value_type(int_test_value))
        @test S isa Pseudoscalar{Float64}
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        S = Pseudoscalar(test_dim, value_type(abs(int_test_value)))
        @test S isa Pseudoscalar{Float64}
    end

    # Bool
    S = Pseudoscalar(test_dim, true)
    @test S isa Pseudoscalar{Float64}

    S = Pseudoscalar(test_dim, false)
    @test S isa Pseudoscalar{Float64}

    # --- Pseudoscalar{T}(dim::Integer, value::Integer) where {T<:AbstractFloat}

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        for precision_type in subtypes(AbstractFloat)
            S = Pseudoscalar{precision_type}(test_dim,
                                             value_type(int_test_value))
            @test S isa Pseudoscalar{precision_type}
        end
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        for precision_type in subtypes(AbstractFloat)
            S = Pseudoscalar{precision_type}(test_dim,
                                             value_type(abs(int_test_value)))
            @test S isa Pseudoscalar{precision_type}
        end
    end

    # Bool
    for precision_type in subtypes(AbstractFloat)
        S = Pseudoscalar{precision_type}(test_dim, true)
        @test S isa Pseudoscalar{precision_type}
    end
end


# --- Function tests

@testset "Pseudoscalar: AbstractBlade interface tests" begin
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

        # value = 0
        S = Pseudoscalar(test_dim, precision_type(0))
        @test dim(S) == test_dim
        @test grade(S) == test_dim
        @test basis(S) == LinearAlgebra.I
        @test volume(S) isa precision_type
        @test volume(S) == 0
        @test norm(S) isa precision_type
        @test norm(S) == 0
        @test sign(S) == 0

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

@testset "Pseudoscalar: Pseudoscalar interface tests" begin
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

@testset "Pseudoscalar: convert() tests" begin
    # Preparations
    test_dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            converted_test_value = precision_type_src(test_value)
            S = Pseudoscalar{precision_type_src}(test_dim, converted_test_value)

            @test convert(Pseudoscalar{precision_type_converted}, S) isa
                  Pseudoscalar{precision_type_converted}
        end
    end
end
