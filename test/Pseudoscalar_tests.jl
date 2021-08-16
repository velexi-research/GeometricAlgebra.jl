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

# --- File inclusions

# Test utilities
include("test_utils.jl")

# --- Constructor tests

@testset "Pseudoscalar: inner constructor" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # --- Pseudoscalar{T}(dim::Integer, value::AbstractFloat)

    test_dim = 10

    # value != 0
    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B isa Pseudoscalar{precision_type}
            @test B.dim == test_dim
            @test B.value isa precision_type
            @test B.value == precision_type(converted_test_value)
        end
    end
    
    # value == 0
    test_value = 0
    
    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            converted_test_value = value_type(test_value)

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B isa Zero{precision_type}
        end
    end

    # --- Pseudoscalar{T}(dim::Integer, value::Integer)

    test_dim = 10

    # value != 0
    test_value = (rand() > 0.5) ? 10 : -10

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            converted_test_value = value_type(test_value)

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B.dim == test_dim
            @test B.value isa precision_type
            @test B.value == precision_type(converted_test_value)
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = value_type(abs(test_value))

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B.dim == test_dim
            @test B.value isa precision_type
            @test B.value == precision_type(converted_test_value)
        end

        # Bool
        B = Pseudoscalar{precision_type}(test_dim, true)
        @test B isa Pseudoscalar{precision_type}
        @test B.dim == test_dim
        @test B.value isa precision_type
        @test B.value == precision_type(1)

        B = Pseudoscalar{precision_type}(test_dim, false)
        @test B isa Zero{precision_type}
    end

    # value == 0
    test_value = 0

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            converted_test_value = value_type(test_value)

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B isa Zero{precision_type}
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            converted_test_value = value_type(abs(test_value))

            B = Pseudoscalar{precision_type}(test_dim, converted_test_value)
            @test B isa Zero{precision_type}
        end
    end

    # --- Invalid data fields

    for precision_type in subtypes(AbstractFloat)
        # dim == 0
        @test_throws ArgumentError Pseudoscalar{precision_type}(0, 10)

        # dim < 0
        @test_throws ArgumentError Pseudoscalar{precision_type}(-1, 10)
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

    # value != 0
    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        converted_test_value = precision_type(test_value)
        B = Pseudoscalar(test_dim, converted_test_value)
        @test B isa Pseudoscalar{precision_type}
        @test B.dim == test_dim
        @test B.value == converted_test_value
    end

    # value == 0
    test_value = 0

    for precision_type in subtypes(AbstractFloat)
        converted_test_value = precision_type(test_value)
        B = Pseudoscalar(test_dim, converted_test_value)
        @test B isa Zero{precision_type}
    end

    # --- Pseudoscalar(dim::Integer, value::Integer)

    # value != 0
    test_value = (rand() > 0.5) ? 10 : -10

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        B = Pseudoscalar(test_dim, value_type(test_value))
        @test B isa Pseudoscalar{Float64}
        @test B.dim == test_dim
        @test B.value == value_type(test_value)
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        B = Pseudoscalar(test_dim, value_type(abs(test_value)))
        @test B isa Pseudoscalar{Float64}
        @test B.dim == test_dim
        @test B.value == value_type(abs(test_value))
    end

    # Bool
    B = Pseudoscalar(test_dim, true)
    @test B isa Pseudoscalar{Float64}
    @test B.dim == test_dim
    @test B.value == 1.0

    B = Pseudoscalar(test_dim, false)
    @test B isa Zero{Float64}
    
    # value == 0
    test_value = 0

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        B = Pseudoscalar(test_dim, value_type(test_value))
        @test B isa Zero{Float64}
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        B = Pseudoscalar(test_dim, value_type(abs(test_value)))
        @test B isa Zero{Float64}
    end

    # --- Invalid data fields

    # dim == 0
    @test_throws ArgumentError Pseudoscalar(0, 10)

    # dim < 0
    @test_throws ArgumentError Pseudoscalar(-1, 10)
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

    test_value = get_random_value()

    # --- Pseudoscalar(dim::Integer, value::T) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)
        B = Pseudoscalar(test_dim, converted_test_value)

        # Construct a Pseudoscalar representing the same pseudoscalar as `B`
        B_copy = Pseudoscalar(B)
        @test B_copy isa Pseudoscalar{precision_type}
        @test B.dim == test_dim
        @test B.value == precision_type(converted_test_value)

        # Construct a Pseudoscalar representing the same space as `B` with a
        # different value.
        B_copy = Pseudoscalar(B, value=converted_test_value + 1)
        @test B_copy isa Pseudoscalar{precision_type}
        @test B.dim == test_dim
        @test value(B_copy) == converted_test_value + 1
    end
end

# --- Test attribute methods

@testset "Pseudoscalar: AbstractMultivector attribute functions" begin
    # Preparations
    test_dim = 10

    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        B = Pseudoscalar(test_dim, positive_test_value)
        @test dim(B) == test_dim
        @test grades(B) == [test_dim]
        @test blades(B) == [B]
        @test norm(B) isa precision_type
        @test norm(B) == positive_test_value

        @test B[0] == []
        @test B[1] == []
        @test B[test_dim] == [B]

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        B = Pseudoscalar(test_dim, negative_test_value)
        @test norm(B) isa precision_type
        @test norm(B) == abs(negative_test_value)

        # value = Inf
        B = Pseudoscalar(test_dim, precision_type(Inf))
        @test norm(B) isa precision_type
        @test norm(B) == Inf

        # value = -Inf
        B = Pseudoscalar(test_dim, precision_type(-Inf))
        @test norm(B) isa precision_type
        @test norm(B) == Inf
    end
end

@testset "Pseudoscalar: AbstractBlade attribute functions" begin
    # Preparations
    test_dim = 10

    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        B = Pseudoscalar(test_dim, positive_test_value)
        @test grade(B) == test_dim
        @test basis(B) == LinearAlgebra.I
        @test volume(B) isa precision_type
        @test volume(B) == positive_test_value
        @test sign(B) == 1

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        B = Pseudoscalar(test_dim, negative_test_value)
        @test volume(B) isa precision_type
        @test volume(B) == negative_test_value
        @test sign(B) == -1

        # value = Inf
        B = Pseudoscalar(test_dim, precision_type(Inf))
        @test volume(B) isa precision_type
        @test volume(B) == precision_type(Inf)
        @test sign(B) == 1

        # value = -Inf
        B = Pseudoscalar(test_dim, precision_type(-Inf))
        @test volume(B) isa precision_type
        @test volume(B) == precision_type(-Inf)
        @test sign(B) == -1
    end
end

@testset "Pseudoscalar: Pseudoscalar attribute functions" begin
    # Preparations
    test_dim = 10

    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_test_value = precision_type(test_value)

        # value > 0
        positive_test_value = converted_test_value > 0 ?
            converted_test_value : -converted_test_value
        B = Pseudoscalar(test_dim, positive_test_value)
        @test value(B) isa precision_type
        @test value(B) == positive_test_value

        # value < 0
        negative_test_value = converted_test_value > 0 ?
            -converted_test_value : converted_test_value
        B = Pseudoscalar(test_dim, negative_test_value)
        @test value(B) isa precision_type
        @test value(B) == negative_test_value

        # value = Inf
        B = Pseudoscalar(test_dim, precision_type(Inf))
        @test value(B) isa precision_type
        @test value(B) == precision_type(Inf)

        # value = -Inf
        B = Pseudoscalar(test_dim, precision_type(-Inf))
        @test value(B) isa precision_type
        @test value(B) == precision_type(-Inf)
    end
end

# --- Test comparison operations

@testset "Pseudoscalar: ==(B, C)" begin
    # Preparations
    dim = 10

    test_value = get_random_value()

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

@testset "Pseudoscalar: ≈(B, C)" begin
    # Preparations
    dim = 10

    test_value = get_random_value()

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

# --- Tests for AbstractMultivector interface functions

@testset "Pseudoscalar: -(B)" begin
    # Preparations
    test_dim = 10

    test_value = get_random_value(2)  # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)

        minus_B = -B
        @test minus_B isa Pseudoscalar{precision_type}
        @test minus_B == Pseudoscalar{precision_type}(test_dim, -test_value)
    end
end

@testset "Pseudoscalar: reverse(B)" begin
    # --- Preparations

    test_value = get_random_value(1)  # add 1 to avoid 0

    for precision_type in subtypes(AbstractFloat)
        for test_dim in 5:8
            B = Pseudoscalar(test_dim, test_value)
            reverse_B = reverse(B)
            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar(test_dim, test_value) :
                Pseudoscalar(test_dim, -test_value)
            @test reverse(B) == expected_result

            @test B * reverse(B) ≈ norm(B)^2
        end
    end
end

@testset "Pseudoscalar: dual(B)" begin
    # Preparations
    test_value = get_random_value(2)  # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        for test_dim in 5:8
            B = Pseudoscalar{precision_type}(test_dim, test_value)

            dual_B = dual(B)
            expected_result = Scalar{precision_type}(precision_type(test_value))
            @test dual_B isa Scalar{precision_type}
            @test dual_B == expected_result

            # Check dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
            if mod(dim(B), 4) < 2
                @test dual(dual_B, dim(B)) == B
            else
                @test dual(dual_B, dim(B)) == -B
            end
        end
    end
end

@testset "Pseudoscalar: convert(B)" begin
    # Preparations
    test_dim = 10

    test_value = get_random_value()

    # Tests
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            converted_test_value = precision_type_src(test_value)
            B = Pseudoscalar{precision_type_src}(test_dim, converted_test_value)

            # Exercise functionality and check results
            B_converted = convert(precision_type_converted, B)

            @test B_converted isa Pseudoscalar{precision_type_converted}

            if precision_type_src == precision_type_converted
                @test B_converted === B
            else
                @test B_converted !== B
                @test B_converted ≈ B
            end
        end
    end
end

# --- Tests for AbstractBlade interface functions

@testset "Pseudoscalar: reciprocal(B)" begin
    # Preparations
    test_value = get_random_value(2)  # add 2 to avoid 0 and 1

    # Tests
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)

        for test_dim = 5:8
            B = Pseudoscalar(test_dim, converted_value)

            reciprocal_B = reciprocal(B)
            @test reciprocal_B isa Pseudoscalar{precision_type}

            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1 / converted_value) :
                Pseudoscalar{precision_type}(test_dim, -1 / converted_value)
            @test reciprocal_B == expected_result

            @test B * reciprocal_B ≈ 1
        end
    end
end
