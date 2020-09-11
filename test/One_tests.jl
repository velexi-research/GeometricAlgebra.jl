"""
Unit tests for the One type.

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
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra


# ---  Tests

@testset "One type: constructor tests" begin
    # One()
    @test One() === One{Float64}()

    # One(B::AbstractBlade{T}) where {T<:AbstractFloat}
    for type in subtypes(AbstractFloat)
        @test One(Blade{type}([1 2 3])) === One{type}()
        @test One(Scalar{type}(1)) === One{type}()
        @test One(One{type}()) === One{type}()
        @test One(One{type}()) === One{type}()
    end

    # One(::Type{T}) where {T<:AbstractFloat}
    for type in subtypes(AbstractFloat)
        @test One(type) === One{type}()
    end

    # One(::Type{T}) where {T<:AbstractBlade}
    for blade_type in (Blade, Scalar, Zero, One)
        @test One(blade_type) === One{Float64}()
    end

    # One(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}}
    for blade_type in (Blade, Scalar, Zero, One)
        for type in subtypes(AbstractFloat)
            @test One(blade_type{type}) === One{type}()
        end
    end
end

@testset "One type: function tests" begin
    # dim()
    for type in subtypes(AbstractFloat)
        @test dim(One(type)) == 0
    end

    # grade()
    for type in subtypes(AbstractFloat)
        @test grade(One(type)) == 0
    end

    # norm()
    for type in subtypes(AbstractFloat)
        @test norm(One(type)) == 1
    end

    # basis()
    for type in subtypes(AbstractFloat)
        @test basis(One(type)) === nothing
    end

    # inverse()
    for type in subtypes(AbstractFloat)
        @test inverse(One(type)) === One(type)
    end
end

@testset "One type: comparison operation tests" begin
    # :(==)
    for type in subtypes(AbstractFloat)
        @test One(type) == 1
        @test 1 == One(type)
    end
    for type1 in subtypes(AbstractFloat)
        for type2 in subtypes(AbstractFloat)
            @test One(type1) == One(type2)
        end
    end
end
