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

@testset "One: constructor tests" begin
    # --- Inner constructor

    for precision_type in subtypes(AbstractFloat)
        S = One{precision_type}()
        @test S isa One{precision_type}
        @test S === One{precision_type}()
    end

    # --- Outer constructor

    S = One()
    @test S isa One{Float64}
end

@testset "one()" begin
    # one(M::AbstractMultivector)
    for precision_type in subtypes(AbstractFloat)
        # --- M::AbstractScalar

        # M::Zero
        S = one(Zero{precision_type}())
        @test S isa One{precision_type}

        # M::One
        S = one(One{precision_type}())
        @test S isa One{precision_type}

        # M::Scalar
        S = one(Scalar{precision_type}(5))
        @test S isa One{precision_type}

        # --- M::AbstractBlade

        # M::Blade
        S = one(Blade{precision_type}([1 2 3]))
        @test S isa One{precision_type}

        # M::Pseudoscalar
        S = one(Pseudoscalar{precision_type}(10, 5))
        @test S isa One{precision_type}

        # --- M::AbstractMultivector

        # M::Multivector
        S = one(Multivector{precision_type}([Scalar(3)]))
        @test S isa One{precision_type}
    end

    # one(::Type{<:AbstractMultivector{T}})
    for precision_type in subtypes(AbstractFloat)
        # --- Type{<:AbstractScalar{T}}

        S = one(AbstractScalar{precision_type})
        @test S isa One{precision_type}

        S = one(Zero{precision_type})
        @test S isa One{precision_type}

        S = one(One{precision_type})
        @test S isa One{precision_type}

        S = one(Scalar{precision_type})
        @test S isa One{precision_type}

        # --- Type{<:AbstractBlade{T}}

        S = one(AbstractBlade{precision_type})
        @test S isa One{precision_type}

        S = one(Blade{precision_type})
        @test S isa One{precision_type}

        S = one(Pseudoscalar{precision_type})
        @test S isa One{precision_type}

        # --- Type{<:AbstractMultivector{T}}

        S = one(AbstractMultivector{precision_type})
        @test S isa One{precision_type}

        S = one(Multivector{precision_type})
        @test S isa One{precision_type}
    end
end

@testset "One: AbstractMultivector interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = One{precision_type}()
        @test dim(S) == 0
        @test grades(S) == [0]
        @test blades(S) == [S]
        @test S[0] == [S]
        @test S[1] == []
        @test norm(S) isa precision_type
        @test norm(S) == 1
    end
end

@testset "One: AbstractBlade interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = One{precision_type}()
        @test grade(S) == 0
        @test basis(S) == 1
        @test volume(S) isa precision_type
        @test volume(S) == 1
        @test sign(S) == 1
    end
end

@testset "One: AbstractScalar interface functions" begin
    for precision_type in subtypes(AbstractFloat)
        S = One{precision_type}()
        @test value(S) isa precision_type
        @test value(S) == 1
    end
end

@testset "One: convert(S)" begin
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            S = One{precision_type_src}()
            S_converted = convert(AbstractScalar{precision_type_converted}, S)
            @test S_converted isa One{precision_type_converted}
        end
    end
end
