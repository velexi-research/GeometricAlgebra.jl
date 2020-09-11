"""
Unit tests for the Scalar type.

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


# --- Constructor tests

@testset "Scalar constructor tests: typeof(value) = AbstractFloat" begin
    # --- Float64

    # nonzero value
    value = Float64(10.)
    B = Scalar(value)
    @test B.value == value
    @test typeof(B.value) == Float64

    # zero value
    value = Float64(0.)
    B = Scalar(value)
    @test B === Zero

    # --- Float32

    # nonzero value
    value = Float32(10.)
    B = Scalar(value)
    @test B.value == value
    @test typeof(B.value) == Float32

    # zero value
    value = Float32(0.)
    B = Scalar(value)
    @test B === Zero

    # --- Float16

    # nonzero value
    value = Float16(10.)
    B = Scalar(value)
    @test B.value == value
    @test typeof(B.value) == Float16

    # zero value
    value = Float16(0.)
    B = Scalar(value)
    @test B === Zero
end

@testset "Scalar constructor tests: typeof(value) = Integer" begin
    # --- Int64

    # nonzero value
    value = Int64(10)
    B = Scalar(value)
    @test B.value == value
    @test typeof(B.value) == Float64

    # zero value
    value = Int64(0)
    B = Scalar{Float32}(value)
    @test B === Zero

    # --- Int32

    # nonzero value
    value = Int32(10)
    B = Scalar{Float16}(value)
    @test B.value == value
    @test typeof(B.value) == Float16

    # zero value
    value = Int32(0)
    B = Scalar(value)
    @test B === Zero

    # --- Int16

    # nonzero value
    value = Int16(10)
    B = Scalar{Float32}(value)
    @test B.value == value
    @test typeof(B.value) == Float32

    # zero value
    value = Int16(0)
    B = Scalar{Float16}(value)
    @test B === Zero
end

# --- Function tests

@testset "Scalar function tests" begin
    # Nonzero value
    value = 10
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == value
    @test basis(B) === nothing
    @test inverse(B) == Scalar(1 / value)

    # Inf
    value = Inf
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == Inf
    @test basis(B) === nothing
    @test inverse(B) === Zero

    # -Inf
    value = -Inf
    B = Scalar(value)
    @test dim(B) == 0
    @test grade(B) == 0
    @test norm(B) == Inf
    @test basis(B) === nothing
    @test inverse(B) === Zero
end
