"""
Unit tests for the One type.

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

@testset "One: constructor tests" begin
    # One()
    @test One() === One{Float64}()

    # One(B::AbstractBlade{T}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test One(Blade{precision_type}([1 2 3])) === One{precision_type}()
        @test One(Scalar{precision_type}(1)) === One{precision_type}()
        @test One(Zero{precision_type}()) === One{precision_type}()
        @test One(One{precision_type}()) === One{precision_type}()
    end

    # One(::Type{T}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test One(precision_type) === One{precision_type}()
    end

    # One(::Type{<:AbstractBlade})
    for blade_type in (Blade, Scalar, Zero, One)
        @test One(blade_type) === One{Float64}()
    end

    # One(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat}
    for blade_type in (Blade, Scalar, Zero, One)
        for precision_type in subtypes(AbstractFloat)
            @test One(blade_type{precision_type}) === One{precision_type}()
        end
    end
end

@testset "One: AbstractBlade interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = One(precision_type)

        # dim()
        @test dim(B) == 0

        # grade()
        @test grade(B) == 0

        # basis()
        @test basis(B) isa precision_type
        @test basis(B) == 1

        # volume()
        @test volume(B) isa precision_type
        @test volume(B) == 1

        # norm()
        @test norm(B) isa precision_type
        @test norm(B) == 1

        # sign()
        @test sign(B) == 1
    end
end

@testset "One: AbstractScalar interface tests" begin
    for precision_type in subtypes(AbstractFloat)
        # Preparations
        B = One(precision_type)

        # value()
        @test value(B) == 1
    end
end
