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
    # No argument
    @test One() === One{Float64}()

    # Argument: instance::Blade{T}
    @test One(Blade{Float64}([1 2 3])) === One{Float64}()
    @test One(Blade{Float32}([1 2 3])) === One{Float32}()
    @test One(Blade{Float16}([1 2 3])) === One{Float16}()

    # Argument: instance::Scalar{T}
    @test One(Scalar{Float64}(1)) === One{Float64}()
    @test One(Scalar{Float32}(1)) === One{Float32}()
    @test One(Scalar{Float16}(1)) === One{Float16}()

    # Argument: instance::Zero{T}
    @test One(Zero{Float64}()) === One{Float64}()
    @test One(Zero{Float32}()) === One{Float32}()
    @test One(Zero{Float16}()) === One{Float16}()

    # Argument: instance::One{T}
    @test One(One{Float64}()) === One{Float64}()
    @test One(One{Float32}()) === One{Float32}()
    @test One(One{Float16}()) === One{Float16}()

    # Argument: concrete type
    @test One(Blade{Float64}) === One{Float64}()
    @test One(Blade{Float32}) === One{Float32}()
    @test One(Blade{Float16}) === One{Float16}()
    @test One(Scalar{Float64}) === One{Float64}()
    @test One(Scalar{Float32}) === One{Float32}()
    @test One(Scalar{Float16}) === One{Float16}()
    @test One(Zero{Float64}) === One{Float64}()
    @test One(Zero{Float32}) === One{Float32}()
    @test One(Zero{Float16}) === One{Float16}()
    @test One(One{Float64}) === One{Float64}()
    @test One(One{Float32}) === One{Float32}()
    @test One(One{Float16}) === One{Float16}()

    # Argument: parameter type
    @test One(Blade) === One{Float64}()
    @test One(Scalar) === One{Float64}()
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


@testset "One type: comparator tests" begin
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
