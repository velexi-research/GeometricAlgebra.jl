"""
Unit tests for the ScalarInf constants.

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
import GeometricAlgebra.ScalarInfPrecisions


# --- Tests

@testset "ScalarInf: inf() tests" begin
    # inf(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        @test inf(Blade{precision_type}([1 2 3])) ===
            ScalarInfPrecisions[precision_type]
        @test inf(Scalar{precision_type}(1)) ===
            ScalarInfPrecisions[precision_type]
    end

    # inf(::Type{Blade{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test inf(Blade{precision_type}) ===
            ScalarInfPrecisions[precision_type]
    end

    # inf(::Type{Scalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test inf(Scalar{precision_type}) ===
            ScalarInfPrecisions[precision_type]
    end

    @test_throws MethodError inf(AbstractScalar)
end

@testset "ScalarInf: AbstractBlade interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = inf(Blade{precision_type})

        # dim()
        @test dim(B) == 0

        # grade()
        @test grade(B) == 0

        # basis()
        @test basis(B) == 1

        # volume()
        @test volume(B) == Inf

        # norm()
        @test norm(B) == Inf

        # sign()
        @test sign(B) == 1
    end
end

@testset "ScalarInf: AbstractScalar interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = inf(Blade{precision_type})

        # value()
        @test value(B) == Inf
    end
end
