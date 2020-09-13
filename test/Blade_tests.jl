"""
Unit tests for the Blade type.

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

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Constructor tests

@testset "Blade: inner constructor tests" begin
    # Notes
    # -----
    # * Test value of constructed instance

    # --- Blade{T}(vectors::Matrix{T};
    #              atol::Real=blade_zero(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # number of vectors <= dimension of column space
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(vectors)
        @test B.dim == 3
        @test B.grade == 2
        @test size(B.basis) == (3, 2)
        for i in size(B.basis, 2)
            @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
        end
        @test B.norm ≈ 5

        # number of vectors > dimension of column space
        vectors = Matrix{precision_type}([1 2 3; 4 5 6])
        B = Blade{precision_type}(vectors)
        @test B === Zero{precision_type}()

        # vectors are linearly dependent
        vectors = Matrix{precision_type}([1 2 1; 1 2 4; 1 2 9])
        B = Blade{precision_type}(vectors)
        @test B === Zero{precision_type}()

        # norm(blade) < atol
        vectors = Matrix{precision_type}([3 3; 4 4; 0 1])
        B = Blade{precision_type}(vectors, atol=6)
        @test B === Zero{precision_type}()
    end

    # --- Blade{T}(vector::Vector{T};
    #              atol::Real=blade_zero(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        vector = [3, 4, 12]

        # vector is a column vector
        col_vector = Vector{precision_type}(vector)
        B = Blade{precision_type}(col_vector)
        @test B.dim == 3
        @test B.grade == 1
        @test B.basis ≈ reshape(col_vector / 13, length(vector), 1)
        @test LinearAlgebra.norm(B.basis) ≈ 1
        @test B.norm ≈ 13

        # vector is a row vector
        row_vector = reshape(Array{precision_type}(vector), 1, length(vector))
        B = Blade{precision_type}(row_vector)
        @test B.dim == 3
        @test B.grade == 1
        @test B.basis ≈ permutedims(row_vector / 13)
        @test LinearAlgebra.norm(B.basis) ≈ 1
        @test B.norm ≈ 13

        # vector is a zero vector
        zero_vector = Array{precision_type}([0. 0. 0.])
        B = Blade(zero_vector)
        @test B === Zero{precision_type}()

        # norm(vector) < default atol
        small_vector = Vector{precision_type}(vector)
        small_vector /= LinearAlgebra.norm(small_vector)
        small_vector *= blade_atol(precision_type) / 2
        B = Blade{precision_type}(small_vector)
        @test B === Zero{precision_type}()

        # norm(vector) < atol
        col_vector = Vector{precision_type}(vector)
        B = Blade{precision_type}(col_vector,
                                  atol=LinearAlgebra.norm(col_vector) + 1)
        @test B === Zero{precision_type}()
    end
end

@testset "Blade: outer constructor tests" begin
    # Notes
    # -----
    # * Test type of constructed instances. Correct construction of instances
    #   is tested by the inner constructor tests.
    #
    # * Test behavior of `atol` argument.

    # --- Preparations

    vectors = [3 3; -4 -4; 0 1]

    # --- Blade(vectors::Array{T};
    #           atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_vectors = Matrix{precision_type}(vectors)

        # norm(blade) > default atol
        B = Blade(converted_vectors)
        @test B isa Blade{precision_type}

        # norm(blade) < default atol
        small_blade = blade_atol(precision_type) * converted_vectors / 6
        B = Blade(small_blade)
        @test B === Zero{precision_type}()

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B === Zero{precision_type}()
    end

    # --- Blade{T}(vectors::Array{<:AbstractFloat};
    #              atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        for value_type in subtypes(AbstractFloat)
            # Preparations
            converted_vectors = Matrix{value_type}(vectors)

            # norm(blade) > default atol
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) < default atol
            small_blade = blade_atol(precision_type) * converted_vectors / 6
            B = Blade{precision_type}(small_blade)
            @test B === Zero{precision_type}()

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B === Zero{precision_type}()
        end
    end

    # --- Blade(vectors::Array{<:Integer}; atol::Real=blade_atol(Float64))

    # subtypes(Signed)
    for value_type in subtypes(Signed)
        # Preparations
        converted_vectors = Matrix{value_type}(vectors)

        # norm(blade) > default atol
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B === Zero{Float64}()
    end

    # subtypes(Unsigned)
    for value_type in subtypes(Unsigned)
        # Preparations
        converted_vectors = Matrix{value_type}(abs.(vectors))

        # norm(blade) > default atol
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B === Zero{Float64}()
    end

    # Bool: norm(blade) < atol
    converted_vectors = Matrix{Bool}(vectors .!= 0)
    B = Blade(converted_vectors)
    @test B isa Blade{Float64}

    # Bool: norm(blade) < atol
    B = Blade(converted_vectors, atol=2)
    @test B === Zero{Float64}()

    # --- Blade{T}(vectors::Array{<:Integer};
    #              atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # subtypes(Signed)
        for value_type in subtypes(Signed)
            # Preparations
            converted_vectors = Matrix{value_type}(vectors)

            # norm(blade) > default atol
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B === Zero{precision_type}()
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            # Preparations
            converted_vectors = Matrix{value_type}(abs.(vectors))

            # norm(blade) > default atol
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B === Zero{precision_type}()
        end

        # Bool: norm(blade) < atol
        converted_vectors = Matrix{Bool}(vectors .!= 0)
        B = Blade{precision_type}(converted_vectors)
        @test B isa Blade{precision_type}

        # Bool: norm(blade) < atol
        B = Blade{precision_type}(converted_vectors, atol=2)
        @test B === Zero{precision_type}()
    end
end

# --- Function tests

@testset "Blade function tests" begin
    # Construct Blade
    vectors = [3 3; 4 4; 0 1]
    expected_dim, expected_grade = size(vectors)
    B = Blade(vectors)

    @test dim(B) == expected_dim
    @test grade(B) == expected_grade
end
