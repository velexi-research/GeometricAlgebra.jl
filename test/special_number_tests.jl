"""
Unit tests for the special number functions.

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


# --- Special number function tests

@testset "zero() tests" begin
    # zero(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        B = zero(Blade{precision_type}([1 2 3]))
        @test B isa Scalar{precision_type}
        @test B.value == 0

        B = zero(Scalar{precision_type}(5))
        @test B isa Scalar{precision_type}
        @test B.value == 0

        B = zero(Pseudoscalar{precision_type}(10, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 0
    end

    # zero(::Type{Blade{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = zero(Blade{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 0
    end

    # zero(::Type{Scalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = zero(Scalar{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 0
    end

    # zero(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = zero(Pseudoscalar{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 0
    end
end

@testset "one() tests" begin
    # one(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        B = one(Blade{precision_type}([1 2 3]))
        @test B isa Scalar{precision_type}
        @test B.value == 1

        B = one(Scalar{precision_type}(5))
        @test B isa Scalar{precision_type}
        @test B.value == 1

        B = one(Pseudoscalar{precision_type}(100, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 1
    end

    # one(::Type{Blade{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = one(Blade{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 1
    end

    # one(::Type{Scalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = one(Scalar{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 1
    end

    # one(::Type{Pseudoscalar{T}}) where {T<:AbstractFloat}
    for precision_type in subtypes(AbstractFloat)
        B = one(Pseudoscalar{precision_type})
        @test B isa Scalar{precision_type}
        @test B.value == 1
    end
end
