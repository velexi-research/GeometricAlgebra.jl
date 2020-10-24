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
        S = Zero{precision_type}()
        @test S isa Zero{precision_type}
        @test S === Zero{precision_type}()
    end

    # --- Outer constructor

    S = Zero()
    @test S isa Zero{Float64}
end

@testset "zero()" begin
    # zero(M::AbstractMultivector)
    for precision_type in subtypes(AbstractFloat)
        # --- M::AbstractScalar

        # M::Zero
        S = zero(Zero{precision_type}())
        @test S isa Zero{precision_type}

        # M::One
        S = zero(One{precision_type}())
        @test S isa Zero{precision_type}

        # M::Scalar
        S = zero(Scalar{precision_type}(5))
        @test S isa Zero{precision_type}

        # --- M::AbstractBlade

        # M::Blade
        S = zero(Blade{precision_type}([1 2 3]))
        @test S isa Zero{precision_type}

        # M::Pseudoscalar
        S = zero(Pseudoscalar{precision_type}(10, 5))
        @test S isa Zero{precision_type}

        # --- M::AbstractMultivector

        # M::Multivector
        S = zero(Multivector{precision_type}([Scalar(3)]))
        @test S isa Zero{precision_type}
    end

    # zero(::Type{<:AbstractMultivector{T}})
    for precision_type in subtypes(AbstractFloat)
        # --- Type{<:AbstractScalar{T}}

        S = zero(AbstractScalar{precision_type})
        @test S isa Zero{precision_type}

        S = zero(Zero{precision_type})
        @test S isa Zero{precision_type}

        S = zero(One{precision_type})
        @test S isa Zero{precision_type}

        S = zero(Scalar{precision_type})
        @test S isa Zero{precision_type}

        # --- Type{<:AbstractBlade{T}}

        S = zero(AbstractBlade{precision_type})
        @test S isa Zero{precision_type}

        S = zero(Blade{precision_type})
        @test S isa Zero{precision_type}

        S = zero(Pseudoscalar{precision_type})
        @test S isa Zero{precision_type}

        # --- Type{<:AbstractMultivector{T}}

        S = zero(AbstractMultivector{precision_type})
        @test S isa Zero{precision_type}

        S = zero(Multivector{precision_type})
        @test S isa Zero{precision_type}
    end
end

@testset "Zero: AbstractMultivector interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = Zero{precision_type}()
        @test dim(S) == 0
        @test grades(S) == [0]
        @test blades(S) == [S]
        @test S[0] == [S]
        @test S[1] == []
        @test norm(S) isa precision_type
        @test norm(S) == 0
    end
end

@testset "Zero: AbstractBlade interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = Zero{precision_type}()
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == 0
        @test sign(S) == 0
    end
end

@testset "Zero: AbstractScalar interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = Zero{precision_type}()
        @test value(S) isa precision_type
        @test value(S) == 0
    end
end

@testset "Zero: convert(S)" begin
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            S = Zero{precision_type_src}()
            S_converted = convert(AbstractScalar{precision_type_converted}, S)
            @test S_converted isa Zero{precision_type_converted}
        end
    end
end
