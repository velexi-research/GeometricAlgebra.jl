"""
Unit tests for the Zero type.

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


# --- Tests

@testset "Zero type: constructor tests" begin
    # Zero()
    @test Zero() === Zero{Float64}()

    # Zero(B::AbstractBlade{T}) where {T<:AbstractFloat}
    for type in subtypes(AbstractFloat)
        @test Zero(Blade{type}([1 2 3])) === Zero{type}()
        @test Zero(Scalar{type}(1)) === Zero{type}()
        @test Zero(Zero{type}()) === Zero{type}()
        @test Zero(One{type}()) === Zero{type}()
    end

    # Zero(::Type{T}) where {T<:AbstractFloat}
    for type in subtypes(AbstractFloat)
        @test Zero(type) === Zero{type}()
    end

    # Zero(::Type{T}) where {T<:AbstractBlade}
    for blade_type in (Blade, Scalar, Zero, One)
        @test Zero(blade_type) === Zero{Float64}()
    end

    # Zero(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}}
    for blade_type in (Blade, Scalar, Zero, One)
        for type in subtypes(AbstractFloat)
            @test Zero(blade_type{type}) === Zero{type}()
        end
    end
end

@testset "Zero type: function tests" begin
    # dim()
    for type in subtypes(AbstractFloat)
        @test dim(Zero(type)) == 0
    end

    # grade()
    for type in subtypes(AbstractFloat)
        @test grade(Zero(type)) == 0
    end

    # norm()
    for type in subtypes(AbstractFloat)
        @test norm(Zero(type)) == 0
    end

    # basis()
    for type in subtypes(AbstractFloat)
        @test basis(Zero(type)) === nothing
    end

    # inverse()
    for type in subtypes(AbstractFloat)
        @test inverse(Zero(type)) == Scalar(Inf)
    end
end

@testset "Zero type: comparison operation tests" begin
    # :(==)
    for type in subtypes(AbstractFloat)
        @test Zero(type) == 0
        @test 0 == Zero(type)
    end
    for type1 in subtypes(AbstractFloat)
        for type2 in subtypes(AbstractFloat)
            @test Zero(type1) == Zero(type2)
        end
    end
end
