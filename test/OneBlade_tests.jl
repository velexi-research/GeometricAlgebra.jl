"""
Unit tests for the OneBlade type.

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


# --- Tests

@testset "OneBlade: one() tests" begin
    # one(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        @test one(Blade{precision_type}([1 2 3])) === OneBlade()
        @test one(Scalar{precision_type}(1)) === OneBlade()
    end
    @test one(ZeroBlade()) === OneBlade()
    @test one(OneBlade()) === OneBlade()

    # one(::Type{<:AbstractBlade})
    for blade_type in (Blade, AbstractScalar, Scalar, ZeroBlade, OneBlade)
        @test one(blade_type) === OneBlade()
    end
    for precision_type in subtypes(AbstractFloat)
        @test one(Blade{precision_type}) === OneBlade()
        @test one(Scalar{precision_type}) === OneBlade()
    end
end

@testset "OneBlade: AbstractBlade interface tests" begin
    # Preparations
    B = OneBlade()

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

@testset "OneBlade: AbstractScalar interface tests" begin
    # Preparations
    B = OneBlade()

    # value()
    @test value(B) == 1
end
