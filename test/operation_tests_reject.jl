"""
Unit tests for the reject(x, y) function

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
using LinearAlgebra: norm, ⋅
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

# ------ M::Multivector

@testset "reject(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "reject(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "reject(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Zero)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Real)" begin
    @test_skip 1
end

@testset "reject(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "reject(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Zero)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Real)" begin
    @test_skip 1
end

@testset "reject(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "reject(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::Zero)" begin
    @test_skip 1
end

@testset "reject(B::Scalar, C::Real)" begin
    @test_skip 1
end


@testset "reject(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "reject(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::One, C::One)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Zero)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Real)" begin
    @test_skip 1
end

@testset "reject(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "reject(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Zero)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Real)" begin
    @test_skip 1
end

@testset "reject(B::Zero, C::Vector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "reject(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Real, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Real, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Real, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Real, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Real, C::Zero)" begin
    @test_skip 1
end

# ------ B::Vector

@testset "reject(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "reject(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "reject(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "reject(B::Vector, C::Scalar)" begin
    @test_skip 1
end

@testset "reject(B::Vector, C::One)" begin
    @test_skip 1
end

@testset "reject(B::Vector, C::Zero)" begin
    @test_skip 1
end
