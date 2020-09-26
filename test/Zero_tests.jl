"""
Unit tests for the Zero type.

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

@testset "Zero: zero() tests" begin
    # zero(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        @test zero(Blade{precision_type}([1 2 3])) === Zero()
        @test zero(Scalar{precision_type}(1)) === Zero()
    end
    @test zero(Zero()) === Zero()
    @test zero(One()) === Zero()

    # zero(::Type{<:AbstractBlade})
    for blade_type in (Blade, AbstractScalar, Scalar, Zero, One)
        @test zero(blade_type) === Zero()
    end
    for precision_type in subtypes(AbstractFloat)
        @test zero(Blade{precision_type}) === Zero()
        @test zero(Scalar{precision_type}) === Zero()
    end
end

@testset "Zero: AbstractBlade interface tests" begin
    # Preparations
    B = Zero()

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

@testset "Zero: AbstractScalar interface tests" begin
    # Preparations
    B = Zero()

    # value()
    @test value(B) == 0
end
