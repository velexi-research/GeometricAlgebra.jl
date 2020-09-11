"""
Unit tests for the Blade type.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
using Test
import LinearAlgebra

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Unit tests

# ------ Constructor tests

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
    @test B == ZeroBlade(0)

    # single column vector
    vector = [3.; 4.; 12.]
    B = Blade(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float64}(vector/13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vector = [3. 4. 12.]
    B = Blade(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float64}(permutedims(vector/13))
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # vectors are linearly dependent ==> blade is zero
    vectors = [1. 2. 1.; 1. 2. 4; 1. 2. 9]
    B = Blade(vectors)
    @test B == ZeroBlade(0)
end

@testset "Blade constructor tests: typeof(vectors) = Array{AbstractFloat}" begin
    # multiple column vectors. number of vectors <= dimension of column space
    vectors = Array{Float32}([3. 3.; 4. 4; 0. 1.])
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
    @test B == ZeroBlade{Float16}(0)

    # single column vector
    vector = Array{Float32}([3.; 4.; 12.])
    B = Blade(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float32}(vector/13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vector = Array{Float16}([3. 4. 12.])
    B = Blade(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float16}(permutedims(vector/13))
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13
end

@testset "Blade constructor tests: typeof(vectors) = Array{Int}" begin
    # multiple column vectors. number of vectors <= dimension of column space
    vectors = [3 3; 4 4; 0 1]
    B = Blade(vectors)
    @test B.dim == 3
    @test B.grade == 2
    @test size(B.basis) == (3, 2)
    for i in size(B.basis, 2)
        @test LinearAlgebra.norm(B.basis[:, i]) ≈ 1
    end
    @test B.norm ≈ 5

    # multiple column vectors. number of vectors > dimension of column space
    vectors = [1 2 3; 4 5 6]
    B = Blade(vectors)
    @test B == ZeroBlade(0)

    # single column vector
    vector = [3; 4; 12]
    B = Blade{Float32}(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float32}(vector/13)
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13

    # single row vector
    vector = [3 4 12]
    B = Blade{Float16}(vector)
    @test B.dim == 3
    @test B.grade == 1
    @test B.basis ≈ Array{Float16}(permutedims(vector/13))
    @test LinearAlgebra.norm(B.basis) ≈ 1
    @test B.norm ≈ 13
end

# ------ Function tests