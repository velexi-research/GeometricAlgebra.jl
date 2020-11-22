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

@testset "Blade: inner constructors" begin
    #=
      Notes
      -----
      * Test value of constructed instance
    =#

    # --- Preparations

    # Select random volume for tests where `volume` != 0
    test_volume = rand() + 1  # add 1 to avoid 0
    test_volume = rand() > 0.5 ? test_volume : -test_volume

    # --- Default constructor with constraint enforcement
    #
    # Blade{T}(dim::Int, grade::Int, basis::Matrix{T}, volume::Real;
    #          atol::Real=blade_atol(T), enforce_constraints::Bool=true,
    #          copy_basis::Bool=true) where {T<:AbstractFloat}

    # Preparations
    test_dim = 3
    test_grade = 2
    test_volume = 5

    for precision_type in subtypes(AbstractFloat)
        # --- Preparations

        test_basis = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        normalized_basis = copy(test_basis)
        for i in 1:size(normalized_basis, 2)
            normalized_basis[:, i] /= LinearAlgebra.norm(normalized_basis[:, i])
        end

        # --- Valid data fields

        # volume > atol
        # default values for atol, enforce_constraints, copy_basis
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume)

        @test B.dim == test_dim
        @test B.grade == test_grade
        @test B.basis == normalized_basis
        @test B.basis !== normalized_basis
        for i in 1:size(B.basis, 2)
            @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
        end
        @test B.volume == test_volume

        # volume < atol
        # default values for enforce_constraints, copy_basis
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume, atol=abs(test_volume) + 1)
        @test iszero(B)

        # dim == 0
        B = Blade{precision_type}(0, test_grade,
                                  Matrix{precision_type}(rand(0, test_grade)),
                                  test_volume)
        @test iszero(B)

        # grade == 0
        B = Blade{precision_type}(test_dim, 0,
                                  Matrix{precision_type}(rand(test_dim, 0)),
                                  test_volume)
        @test iszero(B)

        # dim == grade, det(basis) > 0
        # default values for atol, enforce_constraints, copy_basis
        pseudoscalar_basis = Matrix{precision_type}([2 1; 4 3])
        for i in 1:size(pseudoscalar_basis, 2)
            pseudoscalar_basis[:, i] /=
                LinearAlgebra.norm(pseudoscalar_basis[:, i])
        end
        B = Blade{precision_type}(2, 2, pseudoscalar_basis, test_volume)

        @test B isa Pseudoscalar{precision_type}
        @test B.dim == 2
        @test B.value == test_volume

        # dim == grade, det(basis) < 0
        # default values for atol, enforce_constraints, copy_basis
        pseudoscalar_basis = Matrix{precision_type}([1 2; 3 4])
        for i in 1:size(pseudoscalar_basis, 2)
            pseudoscalar_basis[:, i] /=
                LinearAlgebra.norm(pseudoscalar_basis[:, i])
        end
        B = Blade{precision_type}(2, 2, pseudoscalar_basis, test_volume)

        @test B isa Pseudoscalar{precision_type}
        @test B.dim == 2
        @test B.value == -test_volume

        # dim == grade, det(basis) == 0
        # default values for atol, enforce_constraints, copy_basis
        pseudoscalar_basis = Matrix{precision_type}([1 4; 2 8])
        for i in 1:size(pseudoscalar_basis, 2)
            pseudoscalar_basis[:, i] /=
                LinearAlgebra.norm(pseudoscalar_basis[:, i])
        end
        B = Blade{precision_type}(2, 2, pseudoscalar_basis, test_volume)
        @test iszero(B)

        # copy_basis = false
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume, copy_basis=false)

        @test B.basis === normalized_basis

        # --- Invalid data fields

        # ------ enforce_constraints = true

        # dim < 0
        @test_throws ArgumentError Blade{precision_type}(-5,
                                                          test_grade,
                                                          test_basis,
                                                          test_volume)

        # grade < 0
        @test_throws ArgumentError Blade{precision_type}(test_dim,
                                                          -1,
                                                          test_basis,
                                                          test_volume)

        # dim != size(basis, 1)
        @test_throws DimensionMismatch Blade{precision_type}(test_dim + 1,
                                                             test_grade,
                                                             test_basis,
                                                             test_volume)

        # grade != size(basis, 2)
        @test_throws DimensionMismatch Blade{precision_type}(test_dim,
                                                             test_grade + 1,
                                                             test_basis,
                                                             test_volume)

        # basis not normalized
        @test_throws ArgumentError Blade{precision_type}(test_dim,
                                                          test_grade,
                                                          test_basis,
                                                          test_volume)

        # ------ enforce_constraints = false

        # dim < 0
        B = Blade{precision_type}(-5, test_grade, test_basis, test_volume,
                                  enforce_constraints=false)
        @test B.dim == -5

        # dim != size(basis, 1)
        # grade != size(basis, 2)
        # basis not normalized
        B = Blade{precision_type}(test_dim + 1, test_grade + 1, test_basis,
                                  test_volume, enforce_constraints=false)
        @test B.dim == test_dim + 1
        @test B.grade == test_grade + 1
        @test B.basis == test_basis
        @test B.basis !== test_basis
    end

    # --- Basic constructor with multiple vectors
    #
    # Blade{T}(vectors::Matrix{<:Real};
    #          volume::Union{Real, Nothing}=nothing,
    #          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    # Preparations
    vectors = Matrix([3 3; -4 -4; 0 1])

    for precision_type in subtypes(AbstractFloat)

        # --- vectors::Matrix{<:AbstractFloat}

        for value_type in subtypes(AbstractFloat)
            # Preparations
            converted_vectors = convert(Matrix{value_type}, vectors)

            # number of vectors <= dimension of column space, volume == nothing
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}
            @test B.dim == 3
            @test B.grade == 2
            @test size(B.basis) == (3, 2)
            for i in 1:size(B.basis, 2)
                @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
            end
            @test B.volume ≈ 5

            # number of vectors <= dimension of column space, volume != nothing
            B = Blade{precision_type}(converted_vectors, volume=test_volume)
            @test B isa Blade{precision_type}
            @test B.dim == 3
            @test B.grade == 2
            @test size(B.basis) == (3, 2)
            for i in 1:size(B.basis, 2)
                @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
            end

            F = LinearAlgebra.qr(converted_vectors)
            signed_norm = prod(LinearAlgebra.diag(F.R))
            @test B.volume == sign(signed_norm) * precision_type(test_volume)

            # number of vectors <= dimension of column space
            # volume < default atol
            B = Blade{precision_type}(converted_vectors,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test iszero(B)

            # one of the dimensions of basis is 0
            @test Blade{precision_type}(rand(0, test_dim)) ==
                zero(Blade{precision_type})
            @test Blade{precision_type}(rand(test_dim, 0)) ==
                zero(Blade{precision_type})

            # number of vectors > dimension of column space
            special_case_vectors = Matrix{value_type}([1 2 3; 4 5 6])
            B = Blade{precision_type}(special_case_vectors)
            @test iszero(B)

            # vectors are linearly dependent, volume == nothing
            special_case_vectors = Matrix{value_type}([1 2 1; 1 2 4; 1 2 9])
            B = Blade{precision_type}(special_case_vectors)
            @test iszero(B)

            # vectors are linearly dependent, volume != nothing
            special_case_vectors = Matrix{value_type}([1 2 1; 1 2 4; 1 2 9])
            B = Blade{precision_type}(special_case_vectors, volume=test_volume)
            @test iszero(B)

            # number of vectors == dimension of column space, volume == nothing
            # det(vectors) > 0
            special_case_vectors = Matrix{value_type}([3 5; 4 12])
            B = Blade{precision_type}(special_case_vectors)
            @test B isa Pseudoscalar{precision_type}
            @test B.dim == 2
            @test B.value ≈ LinearAlgebra.det(special_case_vectors)

            # number of vectors == dimension of column space, volume == nothing
            # det(vectors) < 0
            special_case_vectors = Matrix{value_type}([5 3; 12 4])
            B = Blade{precision_type}(special_case_vectors)
            @test B isa Pseudoscalar{precision_type}
            @test B.dim == 2
            @test B.value ≈ LinearAlgebra.det(special_case_vectors)

            # number of vectors == dimension of column space, volume == nothing
            # det(vectors) == 0
            special_case_vectors = Matrix{value_type}([5 1; 5 1])
            B = Blade{precision_type}(special_case_vectors)
            @test iszero(B)

            # number of vectors == dimension of column space, volume != nothing
            special_case_vectors = Matrix{value_type}([3 5; 4 12])
            B = Blade{precision_type}(special_case_vectors, volume=test_volume)
            @test B isa Pseudoscalar{precision_type}
            @test B.dim == 2
            sign_basis = sign(LinearAlgebra.det(convert(Matrix{precision_type},
                                                        special_case_vectors)))
            @test B.value ≈ sign_basis * precision_type(test_volume)
        end

        # --- vectors::Matrix{<:Integer}
        #
        # Notes
        # -----
        # * We only need to verify the correctness of the type of the result
        #   for a few permutations of the arguments because the computations
        #   for the data fields are independent of the type of the `vectors`
        #   argument. The correctness of the data field computations is
        #   verified in the unit tests for the cases when `vectors` has type
        #   Matrix{<:AbstractFloat}.

        # subtypes(Signed)
        for value_type in subtypes(Signed)
            # Preparations
            converted_vectors = convert(Matrix{value_type}, vectors)

            # number of vectors <= dimension of column space
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # number of vectors <= dimension of column space
            # volume < default atol
            B = Blade{precision_type}(converted_vectors,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # vectors are linearly dependent
            special_case_vectors = Matrix{value_type}([1 2 1; 1 2 4; 1 2 9])
            B = Blade{precision_type}(special_case_vectors)
            @test iszero(B)

            # number of vectors == dimension of column space
            special_case_vectors = Matrix{value_type}([3 5; 4 12])
            B = Blade{precision_type}(special_case_vectors)
            @test B isa Pseudoscalar{precision_type}
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            # Preparations
            converted_vectors = convert(Matrix{value_type}, abs.(vectors))

            # number of vectors <= dimension of column space
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # number of vectors <= dimension of column space
            # volume < default atol
            B = Blade{precision_type}(converted_vectors,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # vectors are linearly dependent
            special_case_vectors = Matrix{value_type}([1 2 1; 1 2 4; 1 2 9])
            B = Blade{precision_type}(special_case_vectors)
            @test iszero(B)

            # number of vectors == dimension of column space
            special_case_vectors = Matrix{value_type}([3 5; 4 12])
            B = Blade{precision_type}(special_case_vectors)
            @test B isa Pseudoscalar{precision_type}
        end

        # --- Bool

        # Preparations
        converted_vectors = Matrix{Bool}(vectors .!= 0)

        # number of vectors <= dimension of column space
        B = Blade{precision_type}(converted_vectors)
        @test B isa Blade{precision_type}

        # number of vectors <= dimension of column space
        # volume < default atol
        B = Blade{precision_type}(converted_vectors,
                                  volume=blade_atol(precision_type) / 2)
        @test iszero(B)

        # vectors are linearly dependent
        special_case_vectors = Matrix{Bool}([1 0 1; 0 1 1; 1 0 1])
        B = Blade{precision_type}(special_case_vectors)
        @test iszero(B)

        # number of vectors == dimension of column space
        special_case_vectors = Matrix{Bool}([1 1; 0 1])
        B = Blade{precision_type}(special_case_vectors)
        @test B isa Pseudoscalar{precision_type}
    end

    # --- Basic constructor with a single vector
    #
    # Blade{T}(vector::Vector{<:Real};
    #          volume::Union{Real, Nothing}=nothing,
    #          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        vector = [3, 4, 12]

        # --- vector::Vector{<:AbstractFloat}

        for value_type in subtypes(AbstractFloat)
            # Preparations
            converted_vector = convert(Vector{value_type}, vector)

            # vector is a column vector, volume == nothing
            B = Blade{precision_type}(converted_vector)
            @test B isa Blade{precision_type}
            @test B.dim == 3
            @test B.grade == 1
            @test B.basis ≈ reshape(converted_vector / 13,
                                    length(converted_vector), 1)
            @test LinearAlgebra.norm(B.basis) ≈ 1
            @test B.volume ≈ 13

            # vector is a column vector, volume != nothing
            B = Blade{precision_type}(converted_vector, volume=test_volume)
            @test B isa Blade{precision_type}
            @test B.dim == 3
            @test B.grade == 1
            @test B.basis ≈ reshape(converted_vector / 13,
                                    length(converted_vector), 1)
            @test LinearAlgebra.norm(B.basis) ≈ 1
            @test B.volume ≈ test_volume

            # vector is a column vector, volume < default atol
            B = Blade{precision_type}(converted_vector,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # vector is a row vector
            row_vector = reshape(converted_vector, 1, length(converted_vector))
            B = Blade{precision_type}(row_vector)
            @test B.dim == 3
            @test B.grade == 1
            @test B.basis ≈ permutedims(row_vector / 13)
            @test LinearAlgebra.norm(B.basis) ≈ 1
            @test B.volume ≈ 13

            # vector is a zero vector, volume == nothing
            zero_vector = Array{value_type}([0., 0., 0.])
            B = Blade{precision_type}(zero_vector)
            @test iszero(B)

            # vector is a zero vector, volume != nothing
            zero_vector = Array{value_type}([0., 0., 0.])
            B = Blade{precision_type}(zero_vector, volume=test_volume)
            @test iszero(B)

            # norm(vector) < default atol
            small_vector = Vector{value_type}(converted_vector)
            small_vector /= LinearAlgebra.norm(small_vector)
            small_vector *= blade_atol(precision_type) / 2
            B = Blade{precision_type}(small_vector)
            @test iszero(B)

            # norm(vector) < atol
            B = Blade{precision_type}(
                converted_vector,
                atol=LinearAlgebra.norm(converted_vector) + 1)
            @test iszero(B)

            # length(vector) == 0
            B = Blade{precision_type}(rand(0))
            @test iszero(B)
        end

        # --- vector::Vector{<:Integer}
        #
        # Notes
        # -----
        # * We only need to verify the correctness of the type of the result
        #   for a few permutations of the arguments because the computations
        #   for the data fields are independent of the type of the `vectors`
        #   argument. The correctness of the data field computations is
        #   verified in the unit tests for the cases when `vector` has type
        #   Vector{<:AbstractFloat}.

        # subtypes(Signed)
        for value_type in subtypes(Signed)
            # Preparations
            converted_vector = convert(Vector{value_type}, vector)

            # vector is a column vector
            B = Blade{precision_type}(converted_vector)
            @test B isa Blade{precision_type}

            # vector is a column vector, volume < default atol
            B = Blade{precision_type}(converted_vector,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # vector is a row vector
            row_vector = reshape(converted_vector, 1, length(converted_vector))
            B = Blade{precision_type}(row_vector)
            @test B isa Blade{precision_type}

            # vector is a zero vector
            zero_vector = Array{value_type}([0, 0, 0])
            B = Blade{precision_type}(zero_vector)
            @test iszero(B)

            # norm(vector) < default atol
            small_vector = Vector{value_type}(converted_vector)
            small_vector /= LinearAlgebra.norm(small_vector)
            small_vector *= blade_atol(precision_type) / 2
            B = Blade{precision_type}(small_vector)
        end

        # subtypes(Unsigned)
        for value_type in subtypes(Unsigned)
            # Preparations
            converted_vector = convert(Vector{value_type}, abs.(vector))

            # vector is a column vector
            B = Blade{precision_type}(converted_vector)
            @test B isa Blade{precision_type}

            # vector is a column vector, volume < default atol
            B = Blade{precision_type}(converted_vector,
                                      volume=blade_atol(precision_type) / 2)
            @test iszero(B)

            # vector is a row vector
            row_vector = reshape(converted_vector, 1, length(converted_vector))
            B = Blade{precision_type}(row_vector)
            @test B isa Blade{precision_type}

            # vector is a zero vector
            zero_vector = Array{value_type}([0, 0, 0])
            B = Blade{precision_type}(zero_vector)
            @test iszero(B)

            # norm(vector) < default atol
            small_vector = Vector{value_type}(converted_vector)
            small_vector /= LinearAlgebra.norm(small_vector)
            small_vector *= blade_atol(precision_type) / 2
            B = Blade{precision_type}(small_vector)
        end

        # --- Bool

        # Preparations
        converted_vector = Vector{Bool}(vector .!= 0)

        # vector is a column vector
        B = Blade{precision_type}(converted_vector)
        @test B isa Blade{precision_type}

        # vector is a column vector, volume < default atol
        B = Blade{precision_type}(converted_vector,
                                  volume=blade_atol(precision_type) / 2)
        @test iszero(B)

        # vector is a row vector
        row_vector = reshape(converted_vector, 1, length(converted_vector))
        B = Blade{precision_type}(row_vector)
        @test B isa Blade{precision_type}

        # vector is a zero vector
        zero_vector = Array{Bool}([0, 0, 0])
        B = Blade{precision_type}(zero_vector)
        @test iszero(B)
    end
end

@testset "Blade: outer constructors - basic constructors" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.

      * Test behavior of keyword arguments: `volume`, `atol`, `copy_basis`.
    =#

    # --- Preparations

    vectors = Matrix([3 3; -4 -4; 0 1])
    one_vector = Vector([3; 4; 0])
    pseudoscalar_vectors = Matrix([3 5; 4 12])

    # Select random volume for tests where `volume` != 0
    test_volume = rand() + 1  # add 1 to avoid 0
    test_volume = rand() > 0.5 ? test_volume : -test_volume

    # --- Blade(vectors::Array{T};
    #           volume::Union{Real, Nothing}=nothing,
    #           atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = convert(Matrix{precision_type}, vectors)

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{precision_type}

        F = LinearAlgebra.qr(converted_vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test volume(B) == sign(signed_norm) * precision_type(test_volume)

        # norm(blade) < default atol
        small_blade = blade_atol(precision_type) * converted_vectors / 6
        B = Blade(small_blade)
        @test iszero(B)

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test iszero(B)

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = convert(Vector{precision_type}, one_vector)

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{precision_type}
        @test volume(B) == precision_type(test_volume)

        # norm(blade) < default atol
        small_blade = blade_atol(precision_type) * converted_one_vector / 6
        B = Blade(small_blade)
        @test iszero(B)

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test iszero(B)

        # --- number of vectors == dimension of column space

        # vectors are linearly independent
        converted_vectors = Matrix{precision_type}(pseudoscalar_vectors)
        B = Blade(converted_vectors)
        @test B isa Pseudoscalar{precision_type}

        # vectors are linearly dependent
        converted_vectors = Matrix{precision_type}([2 3 5; 4 6 10; 1 1 2])
        B = Blade(converted_vectors)
        @test iszero(B)
    end

    # --- Blade(vectors::Array{<:Integer};
    #           volume::Union{Real, Nothing}=nothing,
    #           atol::Real=blade_atol(Float64))

    # subtypes(Signed)
    for vectors_type in subtypes(Signed)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = convert(Matrix{vectors_type}, vectors)

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{Float64}

        F = LinearAlgebra.qr(converted_vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test volume(B) == sign(signed_norm) * Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test iszero(B)

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = convert(Vector{vectors_type}, one_vector)

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{Float64}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test iszero(B)

        # --- number of vectors == dimension of column space

        # vectors are linearly independent
        converted_vectors = Matrix{Float64}(pseudoscalar_vectors)
        B = Blade(converted_vectors)
        @test B isa Pseudoscalar{Float64}

        # vectors are linearly dependent
        converted_vectors = Matrix{Float64}([2 3 5; 4 6 10; 1 1 2])
        B = Blade(converted_vectors)
        @test iszero(B)
    end

    # subtypes(Unsigned)
    for vectors_type in subtypes(Unsigned)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = convert(Matrix{vectors_type}, abs.(vectors))

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{Float64}

        F = LinearAlgebra.qr(converted_vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test volume(B) == sign(signed_norm) * Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test iszero(B)

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = convert(Vector{vectors_type}, one_vector)

        # norm(blade) > default atol, volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{Float64}

        # norm(blade) > default atol, volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test iszero(B)

        # --- number of vectors == dimension of column space

        # vectors are linearly independent
        converted_vectors = Matrix{Float64}(pseudoscalar_vectors)
        B = Blade(converted_vectors)
        @test B isa Pseudoscalar{Float64}

        # vectors are linearly dependent
        converted_vectors = Matrix{Float64}([2 3 5; 4 6 10; 1 1 2])
        B = Blade(converted_vectors)
        @test iszero(B)
    end

    # --- Bool

    # --- vectors isa Matrix

    # Preparations
    converted_vectors = Matrix{Bool}(vectors .!= 0)

    # norm(blade) > default atol, volume == nothing
    B = Blade(converted_vectors)
    @test B isa Blade{Float64}

    # norm(blade) > default atol, volume != nothing
    B = Blade(converted_vectors, volume=test_volume)
    @test B isa Blade{Float64}

    F = LinearAlgebra.qr(converted_vectors)
    signed_norm = prod(LinearAlgebra.diag(F.R))
    @test volume(B) == sign(signed_norm) * Float64(test_volume)

    # norm(blade) < atol
    B = Blade(converted_vectors, atol=6)
    @test iszero(B)

    # --- vectors isa Vector

    # Preparations
    converted_one_vector = Vector{Bool}(one_vector .!= 0)

    # norm(blade) > default atol, volume == nothing
    B = Blade(converted_one_vector)
    @test B isa Blade{Float64}

    # norm(blade) > default atol, volume != nothing
    B = Blade(converted_one_vector, volume=test_volume)
    @test B isa Blade{Float64}
    @test volume(B) == Float64(test_volume)

    # norm(blade) < atol
    B = Blade(converted_one_vector, atol=6)
    @test iszero(B)

    # --- number of vectors == dimension of column space

    # vectors are linearly independent
    converted_vectors = Matrix{Bool}([1 1 0; 0 0 1; 0 1 0])
    B = Blade(converted_vectors)
    @test B isa Pseudoscalar{Float64}
end

@testset "Blade: outer constructors - copy constructors" begin
    #=
      Notes
      -----
      * Test type of constructed instances. Correct construction of instances
        is tested by the inner constructor tests.

      * Test behavior of keyword arguments: `volume`, `atol`, `copy_basis`.
    =#

    # --- Preparations

    vectors = Matrix([3 3; -4 -4; 0 1])
    one_vector = Vector([3; 4; 0])

    # Select random volume for tests where `volume` != 0
    test_volume = rand() + 1  # add 1 to avoid 0
    test_volume = rand() > 0.5 ? test_volume : -test_volume

    # --- Blade{T}(B::Blade{<:AbstractFloat};
    #              volume::Real=volume(B),
    #              atol::Real=blade_atol(T),
    #              copy_basis::Bool=false) where {T<:AbstractFloat}

    for precision_type_src in subtypes(AbstractFloat)
        # Preparations
        vectors = Matrix{precision_type_src}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type_src}(vectors)

        for precision_type_converted in subtypes(AbstractFloat)

            # Construct a Blade representing the same blade as `B`
            B_copy = Blade{precision_type_converted}(B)
            @test B_copy.dim == B.dim
            @test B_copy.grade == B.grade

            if precision_type_src == precision_type_converted
                @test B_copy.basis === B.basis
            else
                @test B_copy.basis ≈ B.basis
                @test B_copy.basis !== B.basis
            end

            @test B_copy.volume ≈ B.volume

            # Construct a Blade representing the same space as `B` with
            # specified volume
            B_copy = Blade{precision_type_converted}(B, volume=test_volume)
            @test B_copy.dim == B.dim
            @test B_copy.grade == B.grade

            if precision_type_src == precision_type_converted
                @test B_copy.basis === B.basis
            else
                @test B_copy.basis ≈ B.basis
                @test B_copy.basis !== B.basis
            end

            @test B_copy.volume == precision_type_converted(test_volume)

            # Construct a Blade representing the same blade as `B` containing
            # a copy of the basis (instead of a reference).
            B_copy = Blade{precision_type_converted}(B, copy_basis=true)
            @test B_copy.dim == B.dim
            @test B_copy.grade == B.grade

            if precision_type_src == precision_type_converted
                @test B_copy.basis == B.basis
            else
                @test B_copy.basis ≈ B.basis
            end
            @test B_copy.basis !== B.basis

            @test B_copy.volume == B.volume
        end
    end

    # --- Blade(B::AbstractBlade{T}; volume=volume(B), copy_basis=false)
    #         where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        converted_vectors = Matrix{precision_type}(vectors)
        B = Blade(converted_vectors)

        # Construct a Blade representing the same blade as `B`
        B_copy = Blade(B)
        @test B_copy isa Blade{precision_type}
        @test B_copy.basis === B.basis

        # Construct a Blade representing the same blade as `B` with specified
        # volume
        B_copy = Blade(B, volume=test_volume)
        @test B_copy isa Blade{precision_type}
        @test volume(B_copy) == precision_type(test_volume)

        # Construct Blade representing the same blade as `B` containing
        # a copy of the basis (instead of a reference).
        B_copy = Blade(B, copy_basis=true)
        @test B_copy isa Blade{precision_type}
        @test B_copy.basis == B.basis
        @test B_copy.basis !== B.basis
    end
end

@testset "Blade: outer constructors - Scalar constructors" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Exercise functionality and check results

    @test Blade(test_value) == Scalar(test_value)
    @test Blade(Scalar(test_value)) == Scalar(test_value)
end

# --- Test attribute methods

@testset "Blade: AbstractMultivector attribute functions" begin
    # --- Preparations

    vectors = [3 3; 4 4; 0 1]
    test_dim, test_grade = size(vectors)

    # --- Test basic functions

    for precision_type in subtypes(AbstractFloat)
        # Blade with sign > 0
        B = Blade{precision_type}(vectors)

        @test dim(B) == test_dim
        @test grades(B) == [test_grade]
        @test blades(B) == [B]

        @test norm(B) isa precision_type
        @test norm(B) ≈ 5

        for k in 0:test_dim
            if k == test_grade
                @test B[k] == [B]
            else
                @test B[k] == []
            end
        end
    end
end

@testset "Blade: AbstractBlade attribute functions" begin
    # --- Preparations

    vectors = [3 3; 4 4; 0 1]
    test_grade = size(vectors, 2)

    # --- Test basic functions

    for precision_type in subtypes(AbstractFloat)
        # Blade with sign > 0
        B = Blade{precision_type}(vectors)

        @test grade(B) == test_grade

        F = LinearAlgebra.qr(vectors)
        @test basis(B) ≈ Matrix(F.Q)

        @test volume(B) isa precision_type
        @test volume(B) ≈ 5

        @test sign(B) == 1

        # Blade with sign < 0
        C = Blade(B, volume=-1)
        @test sign(C) == -1
    end
end

# --- Test comparison operations

@testset "Blade: ==(B, C)" begin
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

@testset "Blade: ≈(B, C)" begin
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

# --- Tests for AbstractMultivector interface functions

@testset "Blade: -(B)" begin
    # Preparations
    vectors = Matrix{Float16}([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)

    # B::Blade
    negative_B = Blade(B, volume=-volume(B))
    @test -B == negative_B

    @test -negative_B == B
end

@testset "Blade: reverse(B)" begin
    # mod(grade, 4) == 1
    vectors = Vector([3; 4; 0; 0; 0])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 2
    vectors = Matrix([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    @test reverse(B) == -B
    @test B * reverse(B) ≈ norm(B)^2

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    @test reverse(B) === B
    @test B * reverse(B) ≈ norm(B)^2
end

@testset "Blade: dual(B)" begin
    # --- Preparations

    # Test values
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- B::Blade

    for dim_B in 10:13
        for grade_B in 2:5
            B = Blade(rand(dim_B, grade_B))

            dual_B = dual(B)

            # --- Check dim, grade, and norm

            @test dim(dual_B) == dim_B
            @test grade(dual_B) == dim_B - grade_B
            @test norm(dual_B) == norm(B)

            # --- Check that B and dual(B) are orthogonal complements

            @test LinearAlgebra.norm(transpose(basis(dual_B)) * basis(B)) <
                10 * eps(Float64)

            # --- Check sign(dual(B))

            # Compute sign of I formed from basis(B) and basis(dual(B))
            sign_Q = sign(LinearAlgebra.det(hcat(basis(B), basis(dual_B))))

            # Compute expected_sign
            expected_sign = sign(B) * sign_Q

            # Account for sign of I^{-1} relative to I
            if mod(dim(B), 4) >= 2
                expected_sign = -expected_sign
            end

            # Account for reversals required to eliminate B
            if mod(grade(B), 4) >= 2
                expected_sign = -expected_sign
            end

            @test sign(dual_B) == expected_sign

            # Check dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
            if mod(dim(B), 4) < 2
                @test dual(dual_B) ≈ B
            else
                @test dual(dual_B) ≈ -B
            end
        end
    end
end

@testset "Blade: convert(B)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]

    # Tests
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            converted_vectors = convert.(precision_type_src, vectors)
            B = Blade{precision_type_src}(converted_vectors)

            # Exercise functionality and check results
            B_converted = convert(Blade{precision_type_converted}, B)

            @test B_converted isa Blade{precision_type_converted}

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

@testset "Blade: reciprocal(B)" begin
    # mod(grade, 4) == 1
    vectors = Vector([3; 4; 0; 0; 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 2
    vectors = Matrix([3 3; 4 4; 0 1; 0 0; 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 3
    vectors = Matrix([3 3 3; 4 4 4; 0 1 0; 0 0 1; 0 0 0])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=-1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1

    # mod(grade, 4) == 0
    vectors = Matrix([3 3 3 3; 4 4 4 4; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    B = Blade(vectors)
    expected_reciprocal = Blade(B, volume=1 / volume(B))
    @test reciprocal(B) ≈ expected_reciprocal
    @test B * reciprocal(B) ≈ 1
end
