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
Unit tests for the dot(x, y) function
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

@testset "dot(M::Multivector, N::Multivector)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "dot(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "dot(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Zero)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Real)" begin
    @test_skip 1
end

@testset "dot(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "dot(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Zero)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Real)" begin
    @test_skip 1
end

@testset "dot(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "dot(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::Zero)" begin
    @test_skip 1
end

@testset "dot(B::Scalar, C::Real)" begin
    @test_skip 1
end


@testset "dot(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "dot(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::One, C::One)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Zero)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Real)" begin
    @test_skip 1
end

@testset "dot(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "dot(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Zero)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Real)" begin
    @test_skip 1
end

@testset "dot(B::Zero, C::Vector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "dot(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Real, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Real, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Real, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Real, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Real, C::Zero)" begin
    @test_skip 1
end

# ------ B::Vector

@testset "dot(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "dot(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "dot(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "dot(B::Vector, C::Scalar)" begin
    @test_skip 1
end

@testset "dot(B::Vector, C::One)" begin
    @test_skip 1
end

@testset "dot(B::Vector, C::Zero)" begin
    @test_skip 1
end
