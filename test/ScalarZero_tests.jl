"""
Unit tests for the ScalarZero constants.

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
import GeometricAlgebra.ScalarZeroPrecisions


# --- Tests

@testset "ScalarZero: zero() tests" begin
    # zero(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        @test zero(Blade{precision_type}([1 2 3])) ===
            ScalarZeroPrecisions[precision_type]
        @test zero(Scalar{precision_type}(1)) ===
            ScalarZeroPrecisions[precision_type]
    end

    # zero(::Type{Blade{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test zero(Blade{precision_type}) ===
            ScalarZeroPrecisions[precision_type]
    end

    # zero(::Type{Scalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test zero(Scalar{precision_type}) ===
            ScalarZeroPrecisions[precision_type]
    end

    @test_throws MethodError zero(AbstractScalar)
end

@testset "ScalarZero: AbstractBlade interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = zero(Blade{precision_type})

        # dim()
        @test dim(B) == 0

        # grade()
        @test grade(B) == 0

        # basis()
        @test basis(B) == 1

        # volume()
        @test volume(B) == 0

        # norm()
        @test norm(B) == 0

        # sign()
        @test sign(B) == 0
    end
end

@testset "ScalarZero: AbstractScalar interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = zero(Blade{precision_type})

        # value()
        @test value(B) == 0
    end
end
