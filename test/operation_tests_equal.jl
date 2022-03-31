#   Copyright (c) 2020-2022 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Unit tests for ==(x, y) and !=(x, y)
"""

# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- File inclusions

# Test utilities
include("test_utils.jl")

# --- Tests

# ------ B::Blade

@testset "==(B::Blade, C::Blade)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) == volume(C)
    # basis(B) == basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            if precision_type1 == precision_type2
                @test B == C
            else
                @test B != C
            end
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B != C
        end
    end

    # volume(B) != volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=2*volume(B))
            @test B != C
        end
    end

    # basis(B) != basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 5; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                      volume=test_volume)
            @test B != C
        end
    end
end

@testset "!=(B::Blade, C::Pseudoscalar)" begin
    test_vectors = [3 3; 4 4; 0 1]

    test_dim = size(test_vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::Scalar)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::One)" begin
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Blade, C::Zero)" begin
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Blade, C::Real)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::Vector)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_vector = [1; 2; 3]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Pseudoscalar

@testset "!=(B::Pseudoscalar, C::Blade)" begin
    test_dim = 3
    test_value = 5

    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::Pseudoscalar, C::Pseudoscalar)" begin
    # Preparations
    dim = 10

    test_value = get_random_value(1)  # add 1 to keep value away from 0

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

@testset "==(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "==(B::Pseudoscalar, C::Real)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "==(B::Pseudoscalar, C::Vector)" begin
    test_value = 5

    test_dim_1 = 3
    test_vector_1 = [3; 4; 1]

    test_dim_2 = 1
    test_vector_2 = [test_value]

    for precision_type in subtypes(AbstractFloat)
        # B != C
        B = Pseudoscalar{precision_type}(test_dim_1, test_value)
        C = Vector{precision_type}(test_vector_1)
        @test B != C

        # B == C
        B = Pseudoscalar{precision_type}(test_dim_2, test_value)
        C = Vector{precision_type}(test_vector_2)
        @test B == C
    end
end

# ------ B::Scalar

@testset "!=(B::Scalar, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::Scalar, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10

    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "==(B::Scalar, C::Scalar)" begin
    # Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
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
end

@testset "==(B::Scalar, C::One)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Scalar, C::Zero)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "==(B::Scalar, C::Real)" begin
    # --- Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    int_value::Int = 3
    float64_or_bigfloat = (Float64, BigFloat)

    # --- Tests

    for precision_type in subtypes(AbstractFloat)
        # B::Scalar, C::AbstractFloat
        B = Scalar(precision_type(test_value))
        for value_type in subtypes(AbstractFloat)
            # ==(B, C)
            if precision_type == value_type
                @test B == value_type(test_value)
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test B == value_type(test_value)
            else
                @test B != value_type(test_value)
            end

            # !=(x,y)
            @test B != value_type(2 * test_value)
        end

        # B::Scalar, C::Integer
        B = Scalar(precision_type(int_value))
        for value_type in subtypes(Signed)
            # ==(B, C)
            @test B == value_type(int_value)

            # !=(x,y)
            @test B != value_type(2 * int_value)
        end

        for value_type in subtypes(Unsigned)
            # ==(B, C)
            @test B == value_type(int_value)

            # !=(x,y)
            @test B != value_type(2 * int_value)
        end

        # B::Scalar, C::Bool
        B = Scalar(precision_type(test_value))
        @test B != true
        @test B != false
    end
end

@testset "==(B::Scalar, C::Vector)" begin
    test_value = 5
    test_vector = [3; 4; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::One

@testset "!=(B::One, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::One, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "==(B::One, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::One, C::One)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = One{precision_type1}()
            C = One{precision_type2}()
            if precision_type1 == precision_type2
                @test B === C
            else
                @test B == C
            end
        end
    end
end

@testset "==(B::One, C::Zero)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "==(B::One, C::Real)" begin
    # B::One, C::Integer
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        # B != C
        C = precision_type(5)
        @test B != C

        # B == C
        C = precision_type(1)
        @test B == C
    end

    # B::One, C::Bool
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        # B != C
        @test B != false

        # B == C
        @test B == true
    end
end

@testset "==(B::One, C::Vector)" begin
    test_vector = [3; 4; 0]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Zero

@testset "!=(B::Zero, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::Zero, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "==(B::Zero, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::Zero, C::One)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Zero, C::Zero)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Zero{precision_type1}()
            C = Zero{precision_type2}()
            if precision_type1 == precision_type2
                @test B === C
            else
                @test B == C
            end
        end
    end
end

@testset "==(B::Zero, C::Real)" begin
    # B::Zero, C::Integer
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        # B != C
        C = precision_type(5)
        @test B != C

        # B == C
        C = precision_type(0)
        @test B == C
    end

    # B::Zero, C::Bool
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        # B != C
        @test B != true

        # B == C
        @test B == false
    end
end

@testset "==(B::Zero, C::Vector)" begin
    test_vector_1 = [3; 4; 0]
    test_vector_2 = [0; 0; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        # B != C
        C = Vector{precision_type}(test_vector_1)
        @test B != C

        # B == C
        C = Vector{precision_type}(test_vector_2)
        @test B == C
    end
end

# ------ B::Real

@testset "!=(B::Real, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::Real, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "==(B::Real, C::Scalar)" begin
    # --- Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    int_value::Int = 3
    float64_or_bigfloat = (Float64, BigFloat)

    # --- Tests

    for precision_type in subtypes(AbstractFloat)
        # B::AbstractFloat, C::Scalar
        C = Scalar(precision_type(test_value))
        for value_type in subtypes(AbstractFloat)
            # ==(B, C)
            if precision_type == value_type
                @test value_type(test_value) == C
            elseif precision_type in float64_or_bigfloat &&
                   value_type in float64_or_bigfloat
                @test value_type(test_value) == C
            else
                @test value_type(test_value) != C
            end

            # !=(x,y)
            @test value_type(2 * test_value) != C
        end

        # B::Integer, C::Scalar
        C = Scalar(precision_type(int_value))
        for value_type in subtypes(Signed)
            # ==(B, C)
            @test value_type(int_value) == C

            # !=(x,y)
            @test value_type(2 * int_value) != C
        end

        for value_type in subtypes(Unsigned)
            # ==(B, C)
            @test value_type(int_value) == C

            # !=(x,y)
            @test value_type(2 * int_value) != C
        end

        # B::Bool, C::Scalar
        C = Scalar(precision_type(test_value))
        @test true != C
        @test false != C
    end
end

@testset "==(B::Real, C::One)" begin
    # B::Integer, C::One
    for precision_type in subtypes(AbstractFloat)
        C = One{precision_type}()

        # B != C
        B = precision_type(5)
        @test B != C

        # B == C
        B = precision_type(1)
        @test B == C
    end

    # B::Bool, C::One
    for precision_type in subtypes(AbstractFloat)
        C = One{precision_type}()
        
        # B != C
        @test false != C

        # B == C
        @test true == C
    end
end

@testset "==(B::Real, C::Zero)" begin
    # B::Integer, C::Zero
    for precision_type in subtypes(AbstractFloat)
        C = Zero{precision_type}()

        # B != C
        B = precision_type(5)
        @test B != C

        # B == C
        B = precision_type(0)
        @test B == C
    end

    # B::Bool, C::Zero
    for precision_type in subtypes(AbstractFloat)
        C = Zero{precision_type}()

        # B != C
        @test true != C

        # B == C
        @test false == C
    end
end

# ------ B::Vector

@testset "!=(B::Vector, C::Blade)" begin
    test_vector = [1; 2; 3]
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "==(B::Vector, C::Pseudoscalar)" begin
    test_value = 5

    test_vector_1 = [1; 2; 3]
    test_dim_1 = 3

    test_vector_2 = [test_value]
    test_dim_2 = 1

    for precision_type in subtypes(AbstractFloat)
        # B != C
        B = Vector{precision_type}(test_vector_1)
        C = Pseudoscalar{precision_type}(test_dim_1, test_value)
        @test B != C

        # B == C
        B = Vector{precision_type}(test_vector_2)
        C = Pseudoscalar{precision_type}(test_dim_2, test_value)
        @test B == C
    end
end

@testset "==(B::Vector, C::Scalar)" begin
    test_vector = [1; 2; 3]
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::Vector, C::One)" begin
    test_vector = [1; 2; 3]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Vector, C::Zero)" begin
    # B != C
    test_vector = [1; 2; 3]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test B != C
    end

    # B == C
    test_vector = [0; 0; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test B == C
    end
end
