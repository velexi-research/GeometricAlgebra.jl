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

@testset "zero()" begin
    # Preparations
    test_dim = 50

    # zero(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        B = zero(Blade{precision_type}([1 2 3]))
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == 3

        B = zero(Scalar{precision_type}(test_dim, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim

        B = zero(Pseudoscalar{precision_type}(test_dim, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim
    end

    # zero(::Type{<:AbstractBlade{T}}, dim)
    for precision_type in subtypes(AbstractFloat)
        B = zero(Blade{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim

        B = zero(Scalar{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim

        B = zero(Pseudoscalar{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim
    end

    # zero(::Type{<:AbstractFloat}, dim)
    for precision_type in subtypes(AbstractFloat)
        B = zero(precision_type, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 0
        @test B.dim == test_dim
    end
end

@testset "one()" begin
    # Preparations
    test_dim = 50

    # one(B::AbstractBlade)
    for precision_type in subtypes(AbstractFloat)
        B = one(Blade{precision_type}([1 2 3]))
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == 3

        B = one(Scalar{precision_type}(test_dim, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim

        B = one(Pseudoscalar{precision_type}(test_dim, 5))
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim
    end

    # one(::Type{<:AbstractBlade{T}}, dim)
    for precision_type in subtypes(AbstractFloat)
        B = one(Blade{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim

        B = one(Scalar{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim

        B = one(Pseudoscalar{precision_type}, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim
    end

    # one(::Type{<:AbstractFloat}, dim)
    for precision_type in subtypes(AbstractFloat)
        B = one(precision_type, test_dim)
        @test B isa Scalar{precision_type}
        @test B.value == 1
        @test B.dim == test_dim
    end
end
