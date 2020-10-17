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

@testset "Blade: inner constructor" begin
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
    #          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        test_dim = 3
        test_grade = 2
        test_volume = 5

        basis = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        normalized_basis = copy(basis)
        for i in size(normalized_basis, 2)
            normalized_basis[:, i] /= LinearAlgebra.norm(basis[:, i])
        end

        # valid data fields
        # default values for atol, enforce_constraints, copy_basis
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume)

        @test B.dim == test_dim
        @test B.grade == test_grade
        @test B.basis == normalized_basis
        @test B.basis !== normalized_basis
        for i in size(B.basis, 2)
            @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
        end
        @test B.volume == test_volume

        # valid data fields, volume < atol
        # default values for enforce_constraints, copy_basis
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume, atol=abs(test_volume) + 1)
        @test B == zero(Blade{precision_type})

        # valid data fields, dim == grade
        # default values for atol, enforce_constraints, copy_basis
        pseudoscalar_basis = Matrix{precision_type}([1 2; 3 4])
        for i in size(pseudoscalar_basis, 2)
            pseudoscalar_basis[:, i] /=
                LinearAlgebra.norm(pseudoscalar_basis[:, i])
        end
        B = Blade{precision_type}(2, 2, pseudoscalar_basis, test_volume)

        @test B isa Pseudoscalar{precision_type}
        @test B.dim == 2
        @test B.value == test_volume

        # invalid data fields: dim != size(basis, 1)
        # enforce_constraints = true
        @test_throws DimensionMismatch Blade{precision_type}(test_dim + 1,
                                                             test_grade,
                                                             basis,
                                                             test_volume)

        # invalid data fields: grade != size(basis, 2)
        # enforce_constraints = true
        @test_throws DimensionMismatch Blade{precision_type}(test_dim,
                                                             test_grade + 1,
                                                             basis,
                                                             test_volume)

        # invalid data fields: basis not normalized
        # enforce_constraints = true
        @test_throws ErrorException Blade{precision_type}(test_dim,
                                                          test_grade,
                                                          basis,
                                                          test_volume)

        # invalid data fields
        # enforce_constraints = false
        B = Blade{precision_type}(test_dim + 1, test_grade + 1, basis,
                                  test_volume, enforce_constraints=false)
        @test B.dim == test_dim + 1
        @test B.grade == test_grade + 1
        @test B.basis == basis
        @test B.basis !== basis

        # valid data fields
        # copy_basis = false
        B = Blade{precision_type}(test_dim, test_grade, normalized_basis,
                                  test_volume, copy_basis=false)

        @test B.basis === normalized_basis
    end

    # --- Basic constructor with multiple vectors
    #
    # Blade{T}(vectors::Matrix{T};
    #          volume::Union{Real, Nothing}=nothing,
    #          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # number of vectors <= dimension of column space; volume == nothing
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(vectors)
        @test B.dim == 3
        @test B.grade == 2
        @test size(B.basis) == (3, 2)
        for i in size(B.basis, 2)
            @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
        end
        @test B.volume ≈ 5

        # number of vectors <= dimension of column space; volume != nothing
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(vectors, volume=test_volume)
        @test B.dim == 3
        @test B.grade == 2
        @test size(B.basis) == (3, 2)
        for i in size(B.basis, 2)
            @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
        end

        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test B.volume == sign(signed_norm) * precision_type(test_volume)

        # number of vectors <= dimension of column space; volume < default atol
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(vectors,
                                  volume=blade_atol(precision_type) / 2)
        @test B == zero(Blade{precision_type})

        # number of vectors > dimension of column space
        vectors = Matrix{precision_type}([1 2 3; 4 5 6])
        B = Blade{precision_type}(vectors)
        @test B == zero(Blade{precision_type})

        # vectors are linearly dependent; volume == nothing
        vectors = Matrix{precision_type}([1 2 1; 1 2 4; 1 2 9])
        B = Blade{precision_type}(vectors)
        @test B == zero(Blade{precision_type})

        # vectors are linearly dependent; volume != nothing
        vectors = Matrix{precision_type}([1 2 1; 1 2 4; 1 2 9])
        B = Blade{precision_type}(vectors, volume=test_volume)
        @test B == zero(Blade{precision_type})

        # norm(blade) < atol
        vectors = Matrix{precision_type}([3 3; 4 4; 0 1])
        B = Blade{precision_type}(vectors, atol=6)
        @test B == zero(Blade{precision_type})

        # number of vectors == dimension of column space; volume == nothing
        vectors = Matrix{precision_type}([3 5; 4 12])
        B = Blade{precision_type}(vectors)
        @test B isa Pseudoscalar{precision_type}
        @test B.dim == 2
        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test B.value ≈ signed_norm

        # number of vectors == dimension of column space; volume != nothing
        vectors = Matrix{precision_type}([3 5; 4 12])
        B = Blade{precision_type}(vectors, volume=test_volume)
        @test B isa Pseudoscalar{precision_type}
        @test B.dim == 2
        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test B.value ≈ sign(signed_norm) * precision_type(test_volume)
    end

    # --- Basic constructor with a single vector
    #
    # Blade{T}(vector::Vector{T};
    #          volume::Union{Real, Nothing}=nothing,
    #          atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        vector = [3, 4, 12]

        # vector is a column vector; volume == nothing
        col_vector = Vector{precision_type}(vector)
        B = Blade{precision_type}(col_vector)
        @test B.dim == 3
        @test B.grade == 1
        @test B.basis ≈ reshape(col_vector / 13, length(vector), 1)
        @test LinearAlgebra.norm(B.basis) ≈ 1
        @test B.volume ≈ 13

        # vector is a column vector; volume != nothing
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(col_vector, volume=test_volume)
        @test B.dim == 3
        @test B.grade == 1
        @test B.basis ≈ reshape(col_vector / 13, length(vector), 1)
        @test LinearAlgebra.norm(B.basis) ≈ 1

        F = LinearAlgebra.qr(vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test B.volume == sign(signed_norm) * precision_type(test_volume)

        # vector is a column vector; volume < default atol
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(col_vector,
                                  volume=blade_atol(precision_type) / 2)
        @test B == zero(Blade{precision_type})

        # vector is a row vector
        row_vector = reshape(Array{precision_type}(vector), 1, length(vector))
        B = Blade{precision_type}(row_vector)
        @test B.dim == 3
        @test B.grade == 1
        @test B.basis ≈ permutedims(row_vector / 13)
        @test LinearAlgebra.norm(B.basis) ≈ 1
        @test B.volume ≈ 13

        # vector is a zero vector; volume == nothing
        zero_vector = Array{precision_type}([0., 0., 0.])
        B = Blade(zero_vector)
        @test B == zero(Blade{precision_type})

        # vector is a zero vector; volume != nothing
        zero_vector = Array{precision_type}([0., 0., 0.])
        B = Blade(zero_vector, volume=test_volume)
        @test B == zero(Blade{precision_type})

        # norm(vector) < default atol
        small_vector = Vector{precision_type}(vector)
        small_vector /= LinearAlgebra.norm(small_vector)
        small_vector *= blade_atol(precision_type) / 2
        B = Blade{precision_type}(small_vector)
        @test B == zero(Blade{precision_type})

        # norm(vector) < atol
        col_vector = Vector{precision_type}(vector)
        B = Blade{precision_type}(col_vector,
                                  atol=LinearAlgebra.norm(col_vector) + 1)
        @test B == zero(Blade{precision_type})
    end

    # --- Copy constructor
    #
    # Blade{T}(B::Blade{T};
    #          volume::Real=volume(B),
    #          copy_basis::Bool=false) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        vectors = Matrix{precision_type}([3 -3; 4 -4; 0 1])
        B = Blade{precision_type}(vectors)

        # Construct a Blade representing the same blade as `B`
        B_copy = Blade{precision_type}(B)
        @test B_copy.dim == B.dim
        @test B_copy.grade == B.grade
        @test B_copy.basis === B.basis
        @test B_copy.volume == B.volume

        # Construct a Blade representing the same space as `B` with specified
        # volume
        B_copy = Blade{precision_type}(B, volume=test_volume)
        @test B_copy.dim == B.dim
        @test B_copy.grade == B.grade
        @test B_copy.basis === B.basis
        @test B_copy.volume == precision_type(test_volume)

        # Construct a Blade representing the same blade as `B` containing
        # a copy of the basis (instead of a reference).
        B_copy = Blade{precision_type}(B, copy_basis=true)
        @test B_copy.dim == B.dim
        @test B_copy.grade == B.grade
        @test B_copy.basis == B.basis
        @test B_copy.basis !== B.basis
        @test B_copy.volume == B.volume
    end

    # --- Type conversion constructor
    #
    #   Blade{T}(B::Blade{S};
    #            volume::Real=volume(B),
    #            copy_basis::Bool=false) where {T<:AbstractFloat,
    #                                           S<:AbstractFloat}

    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            vectors = Matrix{precision_type_src}([3 -3; 4 -4; 0 1])
            B = Blade{precision_type_src}(vectors)

            # Convert precision of `B`
            B_converted = Blade{precision_type_converted}(B)
            @test B_converted.dim == B.dim
            @test B_converted.grade == B.grade
            if precision_type_converted == precision_type_src
                @test B_converted.basis === B.basis
            else
                @test B_converted.basis ≈ B.basis
                @test B_converted.basis !== B.basis
            end
            @test B_converted.volume == B.volume

            # Convert precision of `B` with new `volume`
            B_converted = Blade{precision_type_converted}(B, volume=test_volume)
            @test B_converted.dim == B.dim
            @test B_converted.grade == B.grade
            if precision_type_converted == precision_type_src
                @test B_converted.basis === B.basis
            else
                @test B_converted.basis ≈ B.basis
                @test B_converted.basis !== B.basis
            end
            @test B_converted.volume ≈ test_volume

            # Convert precision of `B` using a copy of the basis (instead of
            # a reference).
            B_converted = Blade{precision_type_converted}(B, copy_basis=true)
            @test B_converted.dim == B.dim
            @test B_converted.grade == B.grade
            if precision_type_converted == precision_type_src
                @test B_converted.basis == B.basis
                @test B_converted.basis !== B.basis
            else
                @test B_converted.basis ≈ B.basis
                @test B_converted.basis !== B.basis
            end
            @test B_converted.volume == B.volume
        end
    end
end

@testset "Blade: outer constructor - basic constructors" begin
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
    #           volume::Union{Real, Nothing}=nothing, atol::Real=blade_atol(T))
    #           atol::Real=blade_atol(T)) where {T<:AbstractFloat}
    #         where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = Matrix{precision_type}(vectors)

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{precision_type}

        F = LinearAlgebra.qr(converted_vectors)
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test volume(B) == sign(signed_norm) * precision_type(test_volume)

        # norm(blade) < default atol
        small_blade = blade_atol(precision_type) * converted_vectors / 6
        B = Blade(small_blade)
        @test B == zero(Blade{precision_type})

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B == zero(Blade{precision_type})

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = Vector{precision_type}(one_vector)

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{precision_type}
        @test volume(B) == precision_type(test_volume)

        # norm(blade) < default atol
        small_blade = blade_atol(precision_type) * converted_one_vector / 6
        B = Blade(small_blade)
        @test B == zero(Blade{precision_type})

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test B == zero(Blade{precision_type})

        # --- number of vectors == dimension of column space

        converted_vectors = Matrix{precision_type}(pseudoscalar_vectors)
        B = Blade{precision_type}(converted_vectors)
        @test B isa Pseudoscalar{precision_type}
    end

    # --- Blade{T}(vectors::Array{<:AbstractFloat};
    #              volume::Union{Real, Nothing}=nothing,
    #              atol::Real=blade_atol(T))
    #         where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        for vectors_type in subtypes(AbstractFloat)
            # --- vectors isa Matrix

            # Preparations
            converted_vectors = Matrix{vectors_type}(vectors)

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_vectors, volume=test_volume)
            @test B isa Blade{precision_type}

            F = LinearAlgebra.qr(converted_vectors)
            signed_norm = prod(LinearAlgebra.diag(F.R))
            @test volume(B) == sign(signed_norm) * precision_type(test_volume)

            # norm(blade) < default atol
            small_blade = blade_atol(precision_type) * converted_vectors / 6
            B = Blade{precision_type}(small_blade)
            @test B == zero(Blade{precision_type})

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B == zero(Blade{precision_type})

            # --- vectors isa Vector

            # Preparations
            converted_one_vector = Vector{precision_type}(one_vector)

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_one_vector)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_one_vector, volume=test_volume)
            @test B isa Blade{precision_type}
            @test volume(B) == precision_type(test_volume)

            # norm(blade) < default atol
            small_blade = blade_atol(precision_type) * converted_one_vector / 6
            B = Blade{precision_type}(small_blade)
            @test B == zero(Blade{precision_type})

            # norm(blade) < atol
            B = Blade{precision_type}(converted_one_vector, atol=6)
            @test B == zero(Blade{precision_type})

            # --- number of vectors == dimension of column space

            converted_vectors = Matrix{vectors_type}(pseudoscalar_vectors)
            B = Blade{precision_type}(converted_vectors)
            @test B isa Pseudoscalar{precision_type}
        end
    end

    # --- Blade(vectors::Array{<:Integer};
    #           volume::Union{Real, Nothing}=nothing,
    #           atol::Real=blade_atol(Float64))

    # subtypes(Signed)
    for vectors_type in subtypes(Signed)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = Matrix{vectors_type}(vectors)

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B == zero(Blade{Float64})

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = Vector{vectors_type}(one_vector)

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{Float64}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test B == zero(Blade{Float64})

        # --- number of vectors == dimension of column space

        converted_vectors = Matrix{vectors_type}(pseudoscalar_vectors)
        B = Blade(converted_vectors)
        @test B isa Pseudoscalar{Float64}
    end

    # subtypes(Unsigned)
    for vectors_type in subtypes(Unsigned)
        # --- vectors isa Matrix

        # Preparations
        converted_vectors = Matrix{vectors_type}(abs.(vectors))

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_vectors)
        @test B isa Blade{Float64}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_vectors, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_vectors, atol=6)
        @test B == zero(Blade{Float64})

        # --- vectors isa Vector

        # Preparations
        converted_one_vector = Vector{vectors_type}(abs.(one_vector))

        # norm(blade) > default atol; volume == nothing
        B = Blade(converted_one_vector)
        @test B isa Blade{Float64}

        # norm(blade) > default atol; volume != nothing
        B = Blade(converted_one_vector, volume=test_volume)
        @test B isa Blade{Float64}
        @test volume(B) == Float64(test_volume)

        # norm(blade) < atol
        B = Blade(converted_one_vector, atol=6)
        @test B == zero(Blade{Float64})

        # --- number of vectors == dimension of column space

        converted_vectors = Matrix{vectors_type}(pseudoscalar_vectors)
        B = Blade(converted_vectors)
        @test B isa Pseudoscalar{Float64}
    end

    # --- Bool

    # ------ vectors isa Matrix

    # norm(blade) > default atol; volume == nothing
    converted_vectors = Matrix{Bool}(vectors .!= 0)
    B = Blade(converted_vectors)
    @test B isa Blade{Float64}

    # norm(blade) > default atol; volume != nothing
    B = Blade(converted_vectors, volume=test_volume)
    @test B isa Blade{Float64}
    @test volume(B) == Float64(test_volume)

    # norm(blade) < atol
    B = Blade(converted_vectors, atol=2)
    @test B == zero(Blade{Float64})

    # ------ vectors isa Vector

    # norm(blade) > atol; volume == nothing
    converted_one_vector = Vector{Bool}(one_vector .!= 0)
    B = Blade(converted_one_vector)
    @test B isa Blade{Float64}

    # norm(blade) > default atol; volume != nothing
    converted_one_vector = Vector{Bool}(one_vector .!= 0)
    B = Blade(converted_one_vector, volume=test_volume)
    @test B isa Blade{Float64}
    @test volume(B) == Float64(test_volume)

    # norm(blade) < atol
    B = Blade(converted_one_vector, atol=2)
    @test B == zero(Blade{Float64})

    # --- Blade{T}(vectors::Array{<:Integer};
    #              volume::Union{Real, Nothing}=nothing,
    #              atol::Real=blade_atol(T)) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # --- vectors isa Matrix

        # subtypes(Signed)
        for vectors_type in subtypes(Signed)
            # --- vectors isa Matrix

            # Preparations
            converted_vectors = Matrix{vectors_type}(vectors)

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_vectors, volume=test_volume)
            @test B isa Blade{precision_type}

            F = LinearAlgebra.qr(converted_vectors)
            signed_norm = prod(LinearAlgebra.diag(F.R))
            @test volume(B) == sign(signed_norm) * precision_type(test_volume)

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B == zero(Blade{precision_type})

            # --- vectors isa Vector

            # Preparations
            converted_one_vector = Vector{vectors_type}(one_vector)

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_one_vector)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_one_vector, volume=test_volume)
            @test B isa Blade{precision_type}
            @test volume(B) == precision_type(test_volume)

            # norm(blade) < atol
            B = Blade{precision_type}(converted_one_vector, atol=6)
            @test B == zero(Blade{precision_type})

            # --- number of vectors == dimension of column space

            converted_vectors = Matrix{vectors_type}(pseudoscalar_vectors)
            B = Blade{precision_type}(converted_vectors)
            @test B isa Pseudoscalar{precision_type}
        end

        # subtypes(Unsigned)
        for vectors_type in subtypes(Unsigned)
            # --- vectors isa Matrix

            # Preparations
            converted_vectors = Matrix{vectors_type}(abs.(vectors))

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_vectors)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_vectors, volume=test_volume)
            @test B isa Blade{precision_type}
            @test volume(B) == precision_type(test_volume)

            F = LinearAlgebra.qr(converted_vectors)
            signed_norm = prod(LinearAlgebra.diag(F.R))
            @test volume(B) == sign(signed_norm) * precision_type(test_volume)

            # norm(blade) < atol
            B = Blade{precision_type}(converted_vectors, atol=6)
            @test B == zero(Blade{precision_type})

            # --- vectors isa Vector

            # Preparations
            converted_one_vector = Vector{vectors_type}(one_vector)

            # norm(blade) > default atol; volume == nothing
            B = Blade{precision_type}(converted_one_vector)
            @test B isa Blade{precision_type}

            # norm(blade) > default atol; volume != nothing
            B = Blade{precision_type}(converted_one_vector, volume=test_volume)
            @test B isa Blade{precision_type}
            @test volume(B) == precision_type(test_volume)

            # norm(blade) < atol
            B = Blade{precision_type}(converted_one_vector, atol=6)
            @test B == zero(Blade{precision_type})

            # --- number of vectors == dimension of column space

            converted_vectors = Matrix{vectors_type}(pseudoscalar_vectors)
            B = Blade{precision_type}(converted_vectors)
            @test B isa Pseudoscalar{precision_type}
        end

        # --- Bool

        # ------ vectors isa Matrix

        # norm(blade) > atol; volume == nothing
        converted_vectors = Matrix{Bool}(vectors .!= 0)
        B = Blade{precision_type}(converted_vectors)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol; volume != nothing
        B = Blade{precision_type}(converted_vectors, volume=test_volume)
        @test B isa Blade{precision_type}

        F = LinearAlgebra.qr(Matrix{precision_type}(converted_vectors))
        signed_norm = prod(LinearAlgebra.diag(F.R))
        @test volume(B) == sign(signed_norm) * precision_type(test_volume)

        # norm(blade) < atol
        B = Blade{precision_type}(converted_vectors, atol=2)
        @test B == zero(Blade{precision_type})

        # ------ vectors isa Vector

        # norm(blade) > atol; volume == nothing
        converted_one_vector = Vector{Bool}(one_vector .!= 0)
        B = Blade{precision_type}(converted_one_vector)
        @test B isa Blade{precision_type}

        # norm(blade) > default atol; volume != nothing
        B = Blade{precision_type}(converted_one_vector, volume=test_volume)
        @test B isa Blade{precision_type}
        @test volume(B) == precision_type(test_volume)

        # norm(blade) < atol
        B = Blade{precision_type}(converted_one_vector, atol=2)
        @test B == zero(Blade{precision_type})
    end
end

@testset "Blade: outer constructor - copy constructor" begin
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


# --- Function tests

@testset "AbstractBlade interface: B::Blade" begin
    # --- Preparations

    vectors = [3 3; 4 4; 0 1]
    expected_dim, expected_grade = size(vectors)

    # --- Test basic functions

    for precision_type in subtypes(AbstractFloat)
        # Blade with sign > 0
        B = Blade{precision_type}(vectors)

        @test dim(B) == expected_dim
        @test grade(B) == expected_grade

        F = LinearAlgebra.qr(vectors)
        @test basis(B) ≈ Matrix(F.Q)

        @test volume(B) isa precision_type
        @test volume(B) ≈ 5

        @test norm(B) isa precision_type
        @test norm(B) ≈ 5

        @test sign(B) == 1

        # Blade with sign < 0
        C = Blade(B, volume=-1)
        @test sign(C) == -1
    end
end

@testset "convert(B): B::Blade" begin
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
