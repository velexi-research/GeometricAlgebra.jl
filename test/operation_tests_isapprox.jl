#   Copyright 2020 Velexi Corporation
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
Unit tests for isapprox(x, y) and !isapprox(x, y)
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

@testset "isapprox(B::Blade, C::Blade)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)

    # dim(B) == dim(C), grade(B) == grade(C), volume(B) ≈ volume(C)
    # basis(B) ≈ basis(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors))
            @test B ≈ C
        end
    end

    # dim(B) != dim(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3; 4 4; 0 1; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # grade(B) != grade(C)
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 3 3; 4 4 4; 0 1 0; 0 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors2))
            @test B ≉ C
        end
    end

    # volume(B) ≉ volume(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors))
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=2*volume(B))
            @test B ≉ C
        end
    end

    # basis(B) ≉ basis(C)
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        vectors2 = [3 4; 4 10; 0 1]
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors2),
                      volume=test_volume)
            @test B ≉ C
        end
    end

    # B and C have opposite orientations
    test_volume = 5.0
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade(convert(Array{precision_type1}, vectors),
                      volume=test_volume)
            C = Blade(convert(Array{precision_type2}, vectors),
                      volume=-test_volume)
            @test B ≉ C
        end
    end
end

@testset "!isapprox(B::Blade, C::Pseudoscalar)" begin
    test_vectors = [3 3; 4 4; 0 1]

    test_dim = 3
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Scalar)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::One)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = One{precision_type}()
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Zero)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Zero{precision_type}()
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Real)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = precision_type(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Blade, C::Vector)" begin
    # --- Preparations
    
    test_vector = [3; 4; 1]

    # B ≉ C
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end

    test_vectors = [3; 4; 10]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end

    # B ≈ C
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade{precision_type1}(test_vector)
            C = Vector{precision_type2}(test_vector)
            @test B ≈ C
        end
    end
end

# ------ B::Pseudoscalar

@testset "!isapprox(B::Pseudoscalar, C::Blade)" begin
    test_dim = 3
    test_value = 5

    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

@testset "isapprox(B::Pseudoscalar, C::Pseudoscalar)" begin
    # Preparations
    dim = 10

    test_value = get_random_value(1)  # add 1 to keep value away from 0

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

@testset "isapprox(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = One{precision_type}()
        @test B ≉ C
    end
end

@testset "isapprox(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Zero{precision_type}()
        @test B ≉ C
    end
end

@testset "isapprox(B::Pseudoscalar, C::Real)" begin
    test_dim = 10
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = precision_type(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Pseudoscalar, C::Vector)" begin
    test_value = 5

    test_dim_1 = 3
    test_vector_1 = [3; 4; 1]

    test_dim_2 = 1
    test_vector_2 = [test_value]

    for precision_type1 in subtypes(AbstractFloat)
        # B ≉ C
        B = Pseudoscalar{precision_type1}(test_dim_1, test_value)
        C = Vector{precision_type1}(test_vector_1)
        @test B ≉ C

        # B ≈ C
        B = Pseudoscalar{precision_type1}(test_dim_2, test_value)
        for precision_type2 in subtypes(AbstractFloat)
            C = Vector{precision_type2}(test_vector_2)
            @test B ≈ C
        end
    end
end

# ------ B::Scalar

@testset "!isapprox(B::Scalar, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

@testset "isapprox(B::Scalar, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Scalar, C::Scalar)" begin
    # Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # B ≈ C
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            @test B ≈ C

            # B ≉ C
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value + 1))
            @test B ≉ C
        end
    end
end

@testset "isapprox(B::Scalar, C::One)" begin
    for precision_type in subtypes(AbstractFloat)
        C = One{precision_type}()

        # B ≈ C
        test_value = 1 - (√eps(precision_type)) / 2
        B = Scalar{precision_type}(test_value)
        @test B ≈ C

        # B ≉ C
        test_value = 5
        B = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Scalar, C::Zero)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Zero{precision_type}()
        @test B ≉ C

        # B ≈ C is not being tested as isapprox between a Scalar and a
        # Zero can never return true with the default atol and rtol value
    end
end

@testset "isapprox(B::Scalar, C::Real)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        B = Scalar{precision_type1}(test_value)

        # B ≉ C
        C = precision_type1(test_value + 1)
        @test B ≉ C

        # B ≈ C
        for precision_type2 in subtypes(AbstractFloat)
            C = precision_type2(test_value)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Scalar, C::Vector)" begin
    test_value = 5
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end
end

# ------ B::One

@testset "!isapprox(B::One, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

@testset "isapprox(B::One, C::Pseudoscalar)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::One, C::Scalar)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        # B ≈ C
        test_value = 1 - (√eps(precision_type)) / 2
        C = Scalar{precision_type}(test_value)
        @test B ≈ C

        # B ≉ C
        test_value = 5
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::One, C::One)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = One{precision_type1}()
            C = One{precision_type2}()
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::One, C::Zero)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Zero{precision_type}()
        @test B ≉ C
    end
end

@testset "isapprox(B::One, C::Real)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        B = One{precision_type1}()

        # B ≉ C
        C = precision_type1(test_value)
        @test B ≉ C

        # B ≈ C
        for precision_type2 in subtypes(AbstractFloat)
            C = precision_type2(1)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::One, C::Vector)" begin
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end
end

# ------ B::Zero

@testset "!isapprox(B::Zero, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

@testset "isapprox(B::Zero, C::Pseudoscalar)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Zero, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
        
        # B ≈ C is not being tested as isapprox between a Zero and a
        # Scalar can never return true with the default atol and rtol value
    end
end

@testset "isapprox(B::Zero, C::One)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = One{precision_type}()
        @test B ≉ C
    end
end

@testset "isapprox(B::Zero, C::Zero)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Zero{precision_type1}()
            C = Zero{precision_type2}()
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Zero, C::Real)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        B = Zero{precision_type1}()

        # B ≉ C
        C = precision_type1(test_value)
        @test B ≉ C

        # B ≈ C
        for precision_type2 in subtypes(AbstractFloat)
            C = precision_type2(0)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Zero, C::Vector)" begin
    # B ≉ C
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end

    test_vector = [0.001; 0.002; 0.003]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test !isapprox(B, C, atol=0.001)
    end

    # B ≈ C
    test_vector = [0; 0; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B ≈ C
    end 

    test_vector = [0.001; 0.002; 0.003]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test isapprox(B, C, atol=0.01)
    end
end

# ------ B::Real

@testset "!isapprox(B::Real, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

@testset "isapprox(B::Real, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Real, C::Scalar)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        B = precision_type1(test_value)

        # B ≉ C
        C = Scalar{precision_type1}(test_value + 1)
        @test B ≉ C

        # B ≈ C
        for precision_type2 in subtypes(AbstractFloat)
            C = Scalar{precision_type2}(test_value)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Real, C::One)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        # B ≉ C
        B = precision_type1(test_value)
        C = One{precision_type1}()
        @test B ≉ C

        # B ≈ C
        B = precision_type1(1)
        for precision_type2 in subtypes(AbstractFloat)
            C = One{precision_type2}()
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Real, C::Zero)" begin
    test_value = 5
    for precision_type1 in subtypes(AbstractFloat)
        # B ≉ C
        B = precision_type1(test_value)
        C = Zero{precision_type1}()
        @test B ≉ C

        # B ≈ C
        B = precision_type1(0)
        for precision_type2 in subtypes(AbstractFloat)
            C = Zero{precision_type2}()
            @test B ≈ C
        end
    end
end

# ------ B::Vector

@testset "isapprox(B::Vector, C::Blade)" begin
    # --- Preparations

    test_vector = [3; 4; 1]

    # B ≉ C
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end

    test_vectors = [3; 4; 10]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end

    # B ≈ C
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Vector{precision_type1}(test_vector)
            C = Blade{precision_type2}(test_vector)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Vector, C::Pseudoscalar)" begin
    test_value = 5

    test_dim_1 = 3
    test_vector_1 = [3; 4; 1]

    test_dim_2 = 1
    test_vector_2 = [test_value]

    for precision_type1 in subtypes(AbstractFloat)
        # B ≉ C
        B = Vector{precision_type1}(test_vector_1)
        C = Pseudoscalar{precision_type1}(test_dim_1, test_value)
        @test B ≉ C

        # B ≈ C
        B = Vector{precision_type1}(test_vector_2)
        for precision_type2 in subtypes(AbstractFloat)
            C = Pseudoscalar{precision_type2}(test_dim_2, test_value)
            @test B ≈ C
        end
    end
end

@testset "isapprox(B::Vector, C::Scalar)" begin
    test_vector = [3; 4; 1]
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "isapprox(B::Vector, C::One)" begin
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = One{precision_type}()
        @test B ≉ C
    end
end

@testset "isapprox(B::Vector, C::Zero)" begin
    # B ≉ C
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test B ≉ C
    end

    test_vector = [0.001; 0.002; 0.003]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test !isapprox(B, C, atol=0.001)
    end

    # B ≈ C
    test_vector = [0; 0; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test B ≈ C
    end 

    test_vector = [0.001; 0.002; 0.003]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test isapprox(B, C, atol=0.01)
    end
end
