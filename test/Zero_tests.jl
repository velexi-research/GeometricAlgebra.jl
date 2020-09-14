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

@testset "Zero: constructor tests" begin
    # Zero()
    @test Zero() === Zero{Float64}()

    # Zero(B::AbstractBlade{T}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = Zero(Blade{precision_type}([1 2 3])) === Zero{precision_type}()
        @test Zero(Blade{precision_type}([1 2 3])) === Zero{precision_type}()
        @test Zero(Scalar{precision_type}(1)) === Zero{precision_type}()
        @test Zero(Zero{precision_type}()) === Zero{precision_type}()
        @test Zero(One{precision_type}()) === Zero{precision_type}()
    end

    # Zero(::Type{T}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        @test Zero(precision_type) === Zero{precision_type}()
    end

    # Zero(::Type{<:AbstractBlade})
    for blade_type in (Blade, Scalar, Zero, One)
        @test Zero(blade_type) === Zero{Float64}()
    end

    # Zero(::Type{<:AbstractBlade{T}}) where {T<:AbstractFloat}
    for blade_type in (Blade, Scalar, Zero, One)
        for precision_type in subtypes(AbstractFloat)
            @test Zero(blade_type{precision_type}) === Zero{precision_type}()
        end
    end
end

@testset "Zero: basic function tests" begin
    # dim()
    for precision_type in subtypes(AbstractFloat)
        @test dim(Zero(precision_type)) == 0
    end

    # grade()
    for precision_type in subtypes(AbstractFloat)
        @test grade(Zero(precision_type)) == 0
    end

    # norm()
    for precision_type in subtypes(AbstractFloat)
        @test norm(Zero(precision_type)) == 0
    end

    # basis()
    for precision_type in subtypes(AbstractFloat)
        @test basis(Zero(precision_type)) === nothing
    end

    # value()
    for precision_type in subtypes(AbstractFloat)
        B = Zero(precision_type)
        @test value(B) isa precision_type
        @test value(B) == precision_type(0)
    end
end
