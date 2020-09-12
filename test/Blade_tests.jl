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
using Test
import LinearAlgebra

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Constructor tests

@testset "Blade constructor tests: typeof(vectors) = Array{Float64}" begin
    # multiple column vectors. number of vectors <= dimension of column space
    vectors = [3. 3.; 4. 4; 0. 1.]
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 2
    @test size(B.basis) == (3, 2)
    for i in size(B.basis, 2)
        @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
    end
    @test B.norm ≈ 5

    # multiple column vectors. number of vectors > dimension of column space
    vectors = [1. 2. 3.; 4. 5. 6.]
    B = Blade(vectors)
    @test B === Zero()

    # single column vector
    vectors = [3.; 4.; 12.]
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float64}(vectors / 13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vectors = [3. 4. 12.]
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float64}(permutedims(vectors / 13))
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # multiple column vectors with norm less than atol
    vectors = [3. 3.; 4. 4; 0. 1.]
    B = Blade(vectors, atol=6)
    @test B === Zero()

    # vectors are linearly dependent ==> blade is zero
    vectors = [1. 2. 1.; 1. 2. 4; 1. 2. 9]
    B = Blade(vectors)
    @test B === Zero()

    # vectors is a single zero vector
    vectors = [0. 0. 0.]
    B = Blade(vectors)
    @test B === Zero()
end

@testset "Blade constructor tests: typeof(vectors) = Array{AbstractFloat}" begin
    # multiple column vectors. number of vectors <= dimension of column space
    vectors = Array{Float32}([3. 3.; 4. 4.; 0. 1.])
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 2
    @test size(B.basis) == (3, 2)
    for i in size(B.basis, 2)
        @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
    end
    @test B.norm ≈ 5

    # multiple column vectors. number of vectors > dimension of column space
    vectors = Array{Float16}([1. 2. 3.; 4. 5. 6.])
    B = Blade(vectors)
    @test B === Zero()

    # single column vector
    vectors = Array{Float32}([3.; 4.; 12.])
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float32}(vectors / 13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vectors = Array{Float16}([3. 4. 12.])
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float16}(permutedims(vectors / 13))
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # multiple column vectors with norm less than atol
    vectors = Array{Float32}([3. 3.; 4. 4; 0. 1.])
    B = Blade(vectors, atol=6)
    @test B === Zero()

    # vectors are linearly dependent ==> blade is zero
    vectors = Array{Float16}([1. 2. 1.; 1. 2. 4; 1. 2. 9])
    B = Blade(vectors)
    @test B === Zero()

    # vectors is a single zero vector
    vectors = Vector{Float32}([0.; 0.; 0.])
    B = Blade(vectors)
    @test B === Zero()
end

@testset "Blade constructor tests: typeof(vectors) = Array{Integer}" begin
    # multiple column vectors. number of vectors <= dimension of column space
    vectors = Array{Int64}([3 3; 4 4; 0 1])
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 2
    @test size(B.basis) == (3, 2)
    for i in size(B.basis, 2)
        @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
    end
    @test B.norm ≈ 5

    # multiple column vectors. number of vectors > dimension of column space
    vectors = Array{Int32}([1 2 3; 4 5 6])
    B = Blade(vectors)
    @test B === Zero()

    # single column vector
    vectors = Array{Int16}([3; 4; 12])
    B = Blade{Float32}(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float32}(vectors / 13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vectors = Vector{Int64}([3; 4; 12])
    B = Blade{Float16}(vectors)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float16}(vectors / 13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # multiple column vectors with norm less than atol
    vectors = Array{Int32}([3 3; 4 4; 0 1])
    B = Blade{Float32}(vectors, atol=6)
    @test B === Zero()

    # vectors are linearly dependent ==> blade is zero
    vectors = Array{Int16}([1 2 1; 1 2 4; 1 2 9])
    B = Blade{Float16}(vectors)
    @test B === Zero()

    # vectors is a single zero vector
    vectors = Vector{Int64}([0; 0; 0])
    B = Blade{Float32}(vectors)
    @test B === Zero
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
