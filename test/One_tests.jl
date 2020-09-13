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

@testset "One type: constructor tests" begin
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

@testset "One type: function tests" begin
    # dim()
    for precision_type in subtypes(AbstractFloat)
        @test dim(One(precision_type)) == 0
    end

    # grade()
    for precision_type in subtypes(AbstractFloat)
        @test grade(One(precision_type)) == 0
    end

    # norm()
    for precision_type in subtypes(AbstractFloat)
        @test norm(One(precision_type)) == 1
    end

    # basis()
    for precision_type in subtypes(AbstractFloat)
        @test basis(One(precision_type)) === nothing
    end

    # inverse()
    for precision_type in subtypes(AbstractFloat)
        @test inverse(One(precision_type)) === One(precision_type)
    end
end

@testset "One type: comparison operation tests" begin
    # :(==)
    for precision_type in subtypes(AbstractFloat)
        @test One(precision_type) == 1
        @test 1 == One(precision_type)
    end
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            @test One(precision_type1) == One(precision_type2)
        end
    end
end
