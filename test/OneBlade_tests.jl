"""
Unit tests for the OneBlade constants.

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
import GeometricAlgebra.OneBladePrecisions


# --- Tests

@testset "OneBlade: one() tests" begin
    # one(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        @test one(Blade{precision_type}([1 2 3])) ===
            OneBladePrecisions[precision_type]
        @test one(Scalar{precision_type}(1)) ===
            OneBladePrecisions[precision_type]
    end

    # one(::Type{Blade{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test one(Blade{precision_type}) ===
            OneBladePrecisions[precision_type]
    end

    # one(::Type{Scalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test one(Scalar{precision_type}) ===
            OneBladePrecisions[precision_type]
    end

    @test_throws MethodError one(AbstractScalar)
end

@testset "OneBlade: AbstractBlade interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = one(Blade{precision_type})

        # dim()
        @test dim(B) == 0

        # grade()
        @test grade(B) == 0

        # basis()
        @test basis(B) == 1

        # volume()
        @test volume(B) == 1

        # norm()
        @test norm(B) == 1

        # sign()
        @test sign(B) == 1
    end
end

@testset "OneBlade: AbstractScalar interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = one(Blade{precision_type})

        # value()
        @test value(B) == 1
    end
end
