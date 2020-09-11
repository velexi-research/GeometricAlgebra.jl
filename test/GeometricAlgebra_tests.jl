"""
Unit tests for the GeometricAlgebra module.

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

# GeometricAlgebra.jl
using GeometricAlgebra
import GeometricAlgebra.zero, GeometricAlgebra.one


# --- Zero type tests

@testset "Zero type tests" begin
    # Property tests
    @test Zero() == 0
    @test Zero{Float32}() == 0
    @test Zero{Float16}() == 0

    # Basic blade functions
    B = Zero()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 0
    @test basis(B) === nothing
    @test inverse(B) == Scalar(Inf)

    # Float32 tests
    B = Zero{Float32}()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 0
    @test basis(B) === nothing
    @test inverse(B) === Scalar(Inf32)
    @test inverse(B) == Scalar(Inf)

    # Float16 tests
    B = Zero{Float16}()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 0
    @test basis(B) === nothing
    @test inverse(B) === Scalar(Inf16)
end

@testset "zero() function tests" begin
    # Argument: instance::Blade{T}
    @test zero(Blade{Float64}([1 2 3])) === Zero{Float64}()
    @test zero(Blade{Float32}([1 2 3])) === Zero{Float32}()
    @test zero(Blade{Float16}([1 2 3])) === Zero{Float16}()

    # Argument: instance::Scalar{T}
    @test zero(Scalar{Float64}(1)) === Zero{Float64}()
    @test zero(Scalar{Float32}(1)) === Zero{Float32}()
    @test zero(Scalar{Float16}(1)) === Zero{Float16}()

    # Argument: instance::Zero{T}
    @test zero(Zero{Float64}()) === Zero{Float64}()
    @test zero(Zero{Float32}()) === Zero{Float32}()
    @test zero(Zero{Float16}()) === Zero{Float16}()

    # Argument: instance::One{T}
    @test zero(One{Float64}()) === Zero{Float64}()
    @test zero(One{Float32}()) === Zero{Float32}()
    @test zero(One{Float16}()) === Zero{Float16}()

    # Argument: concrete type
    @test zero(Blade{Float64}) === Zero{Float64}()
    @test zero(Blade{Float32}) === Zero{Float32}()
    @test zero(Blade{Float16}) === Zero{Float16}()
    @test zero(Scalar{Float64}) === Zero{Float64}()
    @test zero(Scalar{Float32}) === Zero{Float32}()
    @test zero(Scalar{Float16}) === Zero{Float16}()
    @test zero(Zero{Float64}) === Zero{Float64}()
    @test zero(Zero{Float32}) === Zero{Float32}()
    @test zero(Zero{Float16}) === Zero{Float16}()
    @test zero(One{Float64}) === Zero{Float64}()
    @test zero(One{Float32}) === Zero{Float32}()
    @test zero(One{Float16}) === Zero{Float16}()

    # Argument: parameter type
    @test zero(Blade) === Zero{Float64}()
    @test zero(Scalar) === Zero{Float64}()
end


# ---  One type tests

@testset "One type tests" begin
    # Basic blade functions
    B = One()
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == 1
    @test basis(B) === nothing
    @test inverse(B) === One()

    # Parameterization tests
    B = One{Float32}()
    @test inverse(B) !== One()
    @test inverse(B) === One{Float32}()
end

@testset "one() function tests" begin
    # Argument: instance::Blade{T}
    @test one(Blade{Float64}([1 2 3])) === One{Float64}()
    @test one(Blade{Float32}([1 2 3])) === One{Float32}()
    @test one(Blade{Float16}([1 2 3])) === One{Float16}()

    # Argument: instance::Scalar{T}
    @test one(Scalar{Float64}(1)) === One{Float64}()
    @test one(Scalar{Float32}(1)) === One{Float32}()
    @test one(Scalar{Float16}(1)) === One{Float16}()

    # Argument: instance::Zero{T}
    @test one(Zero{Float64}()) === One{Float64}()
    @test one(Zero{Float32}()) === One{Float32}()
    @test one(Zero{Float16}()) === One{Float16}()

    # Argument: instance::One{T}
    @test one(One{Float64}()) === One{Float64}()
    @test one(One{Float32}()) === One{Float32}()
    @test one(One{Float16}()) === One{Float16}()

    # Argument: concrete type
    @test one(Blade{Float64}) === One{Float64}()
    @test one(Blade{Float32}) === One{Float32}()
    @test one(Blade{Float16}) === One{Float16}()
    @test one(Scalar{Float64}) === One{Float64}()
    @test one(Scalar{Float32}) === One{Float32}()
    @test one(Scalar{Float16}) === One{Float16}()
    @test one(Zero{Float64}) === One{Float64}()
    @test one(Zero{Float32}) === One{Float32}()
    @test one(Zero{Float16}) === One{Float16}()
    @test one(One{Float64}) === One{Float64}()
    @test one(One{Float32}) === One{Float32}()
    @test one(One{Float16}) === One{Float16}()

    # Argument: parameter type
    @test one(Blade) === One{Float64}()
    @test one(Scalar) === One{Float64}()
end
