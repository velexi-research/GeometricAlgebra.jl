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

@testset "Zero type constructor tests" begin
    # No argument
    @test Zero() === Zero{Float64}()

    # Argument: instance::Blade{T}
    @test Zero(Blade{Float64}([1 2 3])) === Zero{Float64}()
    @test Zero(Blade{Float32}([1 2 3])) === Zero{Float32}()
    @test Zero(Blade{Float16}([1 2 3])) === Zero{Float16}()

    # Argument: instance::Scalar{T}
    @test Zero(Scalar{Float64}(1)) === Zero{Float64}()
    @test Zero(Scalar{Float32}(1)) === Zero{Float32}()
    @test Zero(Scalar{Float16}(1)) === Zero{Float16}()

    # Argument: instance::Zero{T}
    @test Zero(Zero{Float64}()) === Zero{Float64}()
    @test Zero(Zero{Float32}()) === Zero{Float32}()
    @test Zero(Zero{Float16}()) === Zero{Float16}()

    # Argument: instance::One{T}
    @test Zero(One{Float64}()) === Zero{Float64}()
    @test Zero(One{Float32}()) === Zero{Float32}()
    @test Zero(One{Float16}()) === Zero{Float16}()

    # Argument: concrete type
    @test Zero(Blade{Float64}) === Zero{Float64}()
    @test Zero(Blade{Float32}) === Zero{Float32}()
    @test Zero(Blade{Float16}) === Zero{Float16}()
    @test Zero(Scalar{Float64}) === Zero{Float64}()
    @test Zero(Scalar{Float32}) === Zero{Float32}()
    @test Zero(Scalar{Float16}) === Zero{Float16}()
    @test Zero(Zero{Float64}) === Zero{Float64}()
    @test Zero(Zero{Float32}) === Zero{Float32}()
    @test Zero(Zero{Float16}) === Zero{Float16}()
    @test Zero(One{Float64}) === Zero{Float64}()
    @test Zero(One{Float32}) === Zero{Float32}()
    @test Zero(One{Float16}) === Zero{Float16}()

    # Argument: parameter type
    @test Zero(Blade) === Zero{Float64}()
    @test Zero(Scalar) === Zero{Float64}()
end

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


# ---  One type tests

@testset "One type constructor tests" begin
    # No argument
    @test One() === One{Float64}()

    # Argument: instance::Blade{T}
    @test One(Blade{Float64}([1 2 3])) === One{Float64}()
    @test One(Blade{Float32}([1 2 3])) === One{Float32}()
    @test One(Blade{Float16}([1 2 3])) === One{Float16}()

    # Argument: instance::Scalar{T}
    @test One(Scalar{Float64}(1)) === One{Float64}()
    @test One(Scalar{Float32}(1)) === One{Float32}()
    @test One(Scalar{Float16}(1)) === One{Float16}()

    # Argument: instance::Zero{T}
    @test One(Zero{Float64}()) === One{Float64}()
    @test One(Zero{Float32}()) === One{Float32}()
    @test One(Zero{Float16}()) === One{Float16}()

    # Argument: instance::One{T}
    @test One(One{Float64}()) === One{Float64}()
    @test One(One{Float32}()) === One{Float32}()
    @test One(One{Float16}()) === One{Float16}()

    # Argument: concrete type
    @test One(Blade{Float64}) === One{Float64}()
    @test One(Blade{Float32}) === One{Float32}()
    @test One(Blade{Float16}) === One{Float16}()
    @test One(Scalar{Float64}) === One{Float64}()
    @test One(Scalar{Float32}) === One{Float32}()
    @test One(Scalar{Float16}) === One{Float16}()
    @test One(Zero{Float64}) === One{Float64}()
    @test One(Zero{Float32}) === One{Float32}()
    @test One(Zero{Float16}) === One{Float16}()
    @test One(One{Float64}) === One{Float64}()
    @test One(One{Float32}) === One{Float32}()
    @test One(One{Float16}) === One{Float16}()

    # Argument: parameter type
    @test One(Blade) === One{Float64}()
    @test One(Scalar) === One{Float64}()
end

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
