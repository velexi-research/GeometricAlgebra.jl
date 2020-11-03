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

# --- Constructor tests

@testset "Zero: constructor tests" begin
    # --- Inner constructor

    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test B isa Zero{precision_type}
        @test B === Zero{precision_type}()
    end

    # --- Outer constructor

    B = Zero()
    @test B isa Zero{Float64}
end

@testset "Zero: zero()" begin
    # zero(M::AbstractMultivector)
    for precision_type in subtypes(AbstractFloat)
        # --- M::AbstractScalar

        # M::Zero
        B = zero(Zero{precision_type}())
        @test B isa Zero{precision_type}

        # M::One
        B = zero(One{precision_type}())
        @test B isa Zero{precision_type}

        # M::Scalar
        B = zero(Scalar{precision_type}(5))
        @test B isa Zero{precision_type}

        # --- M::AbstractBlade

        # M::Blade
        B = zero(Blade{precision_type}([1 2 3]))
        @test B isa Zero{precision_type}

        # M::Pseudoscalar
        B = zero(Pseudoscalar{precision_type}(10, 5))
        @test B isa Zero{precision_type}

        # --- M::AbstractMultivector

        # M::Multivector
        B = zero(Multivector{precision_type}([Scalar(3)]))
        @test B isa Zero{precision_type}
    end

    # zero(::Type{<:AbstractMultivector{T}})
    for precision_type in subtypes(AbstractFloat)
        # --- Type{<:AbstractScalar{T}}

        B = zero(AbstractScalar{precision_type})
        @test B isa Zero{precision_type}

        B = zero(Zero{precision_type})
        @test B isa Zero{precision_type}

        B = zero(One{precision_type})
        @test B isa Zero{precision_type}

        B = zero(Scalar{precision_type})
        @test B isa Zero{precision_type}

        # --- Type{<:AbstractBlade{T}}

        B = zero(AbstractBlade{precision_type})
        @test B isa Zero{precision_type}

        B = zero(Blade{precision_type})
        @test B isa Zero{precision_type}

        B = zero(Pseudoscalar{precision_type})
        @test B isa Zero{precision_type}

        # --- Type{<:AbstractMultivector{T}}

        B = zero(AbstractMultivector{precision_type})
        @test B isa Zero{precision_type}

        B = zero(Multivector{precision_type})
        @test B isa Zero{precision_type}
    end
end

# --- Test attribute methods

@testset "Zero: AbstractMultivector attribute functions" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test dim(B) == 0
        @test grades(B) == [0]
        @test blades(B) == [B]
        @test norm(B) isa precision_type
        @test norm(B) == 0

        @test B[0] == [B]
        @test B[1] == []
    end
end

@testset "Zero: AbstractBlade attribute functions" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test grade(B) == 0
        @test basis(B) == 1
        @test volume(B) isa precision_type
        @test volume(B) == 0
        @test sign(B) == 0
    end
end

@testset "Zero: AbstractScalar attribute functions" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test value(B) isa precision_type
        @test value(B) == 0
    end
end

# --- Tests for AbstractMultivector interface functions

@testset "Zero: -(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        @test -B === B
    end
end

@testset "Zero: reverse(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test reverse(B) === B
    end
end

@testset "Zero: dual(B)" begin
    B = Zero()
    expected_message = "The dual of Zero is not well-defined"
    @test_throws ErrorException(expected_message) dual(B)
end

@testset "Zero: convert(B)" begin
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            B = Zero{precision_type_src}()
            B_converted = convert(AbstractScalar{precision_type_converted}, B)
            @test B_converted isa Zero{precision_type_converted}
        end
    end
end

# --- Tests for AbstractBlade interface functions

@testset "Zero: reciprocal(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Scalar{precision_type}
        @test reciprocal_B == Scalar{precision_type}(Inf)
    end
end
