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
    # No argument
    @test Zero() === Zero{Float64}()

    # Argument: instance::Blade{T}
    @test Zero(Blade{Float64}([1 2 3])) === Zero{Float64}()
    @test Zero(Blade{Float32}([1 2 3])) === Zero{Float32}()
    @test Zero(Blade{Float16}([1 2 3])) === Zero{Float16}()

    # Argument: instance::Scalar{T}
    @test Zero(Scalar{Float64}(1)) === Zero{Float64}()
    @test Zero(Scalar{Float32}(1)) === Zero{Float32}()
    @test Zero(Scalar{Float16}(1)) === Zero{Float16}()

    # Argument: instance::Zero{T}
    @test Zero(Zero{Float64}()) === Zero{Float64}()
    @test Zero(Zero{Float32}()) === Zero{Float32}()
    @test Zero(Zero{Float16}()) === Zero{Float16}()

    # Argument: instance::One{T}
    @test Zero(One{Float64}()) === Zero{Float64}()
    @test Zero(One{Float32}()) === Zero{Float32}()
    @test Zero(One{Float16}()) === Zero{Float16}()

    # Argument: concrete type
    @test Zero(Blade{Float64}) === Zero{Float64}()
    @test Zero(Blade{Float32}) === Zero{Float32}()
    @test Zero(Blade{Float16}) === Zero{Float16}()
    @test Zero(Scalar{Float64}) === Zero{Float64}()
    @test Zero(Scalar{Float32}) === Zero{Float32}()
    @test Zero(Scalar{Float16}) === Zero{Float16}()
    @test Zero(Zero{Float64}) === Zero{Float64}()
    @test Zero(Zero{Float32}) === Zero{Float32}()
    @test Zero(Zero{Float16}) === Zero{Float16}()
    @test Zero(One{Float64}) === Zero{Float64}()
    @test Zero(One{Float32}) === Zero{Float32}()
    @test Zero(One{Float16}) === Zero{Float16}()

    # Argument: parameter type
    @test Zero(Blade) === Zero{Float64}()
    @test Zero(Scalar) === Zero{Float64}()
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


@testset "Zero type: comparator tests" begin
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
