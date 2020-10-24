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

@testset "zero()" begin
    # zero(B::AbstractMultivector)
    for precision_type in subtypes(AbstractFloat)
        # --- B::AbstractScalar

        # B::Zero
        B = zero(Zero{precision_type}())
        @test B isa Zero{precision_type}

        # B::One
        B = zero(One{precision_type}())
        @test B isa Zero{precision_type}

        # B::Scalar
        B = zero(Scalar{precision_type}(5))
        @test B isa Zero{precision_type}

        # --- B::AbstractBlade

        # B::Blade
        B = zero(Blade{precision_type}([1 2 3]))
        @test B isa Zero{precision_type}

        # B::Pseudoscalar
        B = zero(Pseudoscalar{precision_type}(10, 5))
        @test B isa Zero{precision_type}

        # --- B::AbstractMultivector

        # B::Multivector
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

@testset "Zero: AbstractMultivector interface functions" begin
    # --- Preparations

    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test dim(B) == 0
        @test grades(B) == [0]
        @test blades(B) == [B]
        @test B[0] == [B]
        @test B[1] == []
        @test norm(B) == 0
        @test norm(B) isa precision_type
    end
end

@testset "Zero: AbstractBlade interface functions" begin
    # --- Preparations

    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test grade(B) == 0
        @test basis(B) == 1
        @test volume(B) == 0
        @test volume(B) isa precision_type
        @test sign(B) == 0
    end
end

@testset "Zero: AbstractScalar interface functions" begin
    # --- Preparations

    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        @test value(B) == 0
        @test value(B) isa precision_type
    end
end
