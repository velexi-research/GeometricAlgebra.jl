"""
Unit tests for !isapprox(x, y)

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
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

# ------ M::Multivector

@testset "!isapprox(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "!isapprox(B::Blade, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "!isapprox(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "!isapprox(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

# ------ B::One

@testset "!isapprox(B::One, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "!isapprox(B::Zero, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "!isapprox(B::Real, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Vector

@testset "!isapprox(B::Vector, M::Multivector)" begin
    @test_skip 1
end
