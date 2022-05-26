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
Unit tests for !=(x, y)
"""

# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Tests

# ------ M::Multivector

@testset "!=(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "!=(B::Blade, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "!=(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "!=(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

# ------ B::One

@testset "!=(B::One, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "!=(B::Zero, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "!=(B::Real, M::Multivector)" begin
    @test_skip 1
end

# ------ B::Vector

@testset "!=(B::Vector, M::Multivector)" begin
    @test_skip 1
end
