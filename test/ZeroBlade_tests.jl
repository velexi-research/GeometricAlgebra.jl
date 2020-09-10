"""
Unit tests for the ZeroBlade type.

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

@testset "ZeroBlade constructor tests: typeof(value) = AbstractFloat" begin
    # --- Float64

    # nonzero value
    value = Float64(10.)
    B = ZeroBlade(value)
    @test B.value == value
    @test typeof(B.value) == Float64

    # zero value
    value = Float64(0.)
    B = ZeroBlade(value)
    @test B.value == 0
    @test typeof(B.value) == Float64

    # --- Float32

    # nonzero value
    value = Float32(10.)
    B = ZeroBlade(value)
    @test B.value == value
    @test typeof(B.value) == Float32

    # zero value
    value = Float32(0.)
    B = ZeroBlade(value)
    @test B.value == 0
    @test typeof(B.value) == Float32

    # --- Float16

    # nonzero value
    value = Float16(10.)
    B = ZeroBlade(value)
    @test B.value == value
    @test typeof(B.value) == Float16

    # zero value
    value = Float16(0.)
    B = ZeroBlade(value)
    @test B.value == 0
    @test typeof(B.value) == Float16
end

@testset "ZeroBlade constructor tests: typeof(value) = Integer" begin
    # --- Int64

    # nonzero value
    value = Int64(10.)
    B = ZeroBlade(value)
    @test B.value == value
    @test typeof(B.value) == Float64

    # zero value
    value = Int64(0.)
    B = ZeroBlade{Float32}(value)
    @test B.value == 0
    @test typeof(B.value) == Float32

    # --- Int32

    # nonzero value
    value = Int32(10.)
    B = ZeroBlade{Float16}(value)
    @test B.value == value
    @test typeof(B.value) == Float16

    # zero value
    value = Int32(0.)
    B = ZeroBlade(value)
    @test B.value == 0
    @test typeof(B.value) == Float64

    # --- Int16

    # nonzero value
    value = Int16(10.)
    B = ZeroBlade{Float32}(value)
    @test B.value == value
    @test typeof(B.value) == Float32

    # zero value
    value = Int16(0.)
    B = ZeroBlade{Float16}(value)
    @test B.value == 0
    @test typeof(B.value) == Float16
end
