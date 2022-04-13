#   Copyright (c) 2020-2022 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Unit tests for the Zero type.
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
    # --- zero(M::AbstractMultivector)

    for precision_type in subtypes(AbstractFloat)
        # --- M::AbstractScalar

        B = zero(Zero{precision_type}())
        @test B isa Zero{precision_type}

        B = zero(One{precision_type}())
        @test B isa Zero{precision_type}

        B = zero(Scalar{precision_type}(5))
        @test B isa Zero{precision_type}

        # --- M::AbstractBlade

        B = zero(Blade{precision_type}([1 2 3]))
        @test B isa Zero{precision_type}

        B = zero(Pseudoscalar{precision_type}(10, 5))
        @test B isa Zero{precision_type}

        # --- M::AbstractMultivector

        B = zero(Multivector{precision_type}([Scalar(3)]))
        @test B isa Zero{precision_type}
    end

    # --- zero(::Type{<:AbstractMultivector})

    B = zero(AbstractScalar)
    @test B isa Zero{Float64}

    B = zero(Zero)
    @test B isa Zero{Float64}

    B = zero(One)
    @test B isa Zero{Float64}

    B = zero(Scalar)
    @test B isa Zero{Float64}

    B = zero(AbstractBlade)
    @test B isa Zero{Float64}

    B = zero(Blade)
    @test B isa Zero{Float64}

    B = zero(Pseudoscalar)
    @test B isa Zero{Float64}

    B = zero(AbstractMultivector)
    @test B isa Zero{Float64}

    B = zero(Multivector)
    @test B isa Zero{Float64}

    # --- zero(::Type{<:AbstractMultivector{T}})
    
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

# --- Test comparison operations

@testset "Zero: iszero(B)" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        # Blade
        vectors = [3 3; 4 4; 0 1; 0 1]
        B = Blade{precision_type}(vectors)
        @test !iszero(B)

        vectors = [1 2; 1 2; 0 0]
        B = Blade{precision_type}(vectors)
        @test iszero(B)

        # Pseudoscalar
        B = Pseudoscalar{precision_type}(10, 3)
        @test !iszero(B)

        # Scalar
        B = Scalar{precision_type}(0)
        @test iszero(B)

        B = Scalar{precision_type}(3)
        @test !iszero(B)

        # Zero
        B = Zero{precision_type}()
        @test iszero(B)

        # One
        B = One{precision_type}()
        @test !iszero(B)
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

@testset "Zero: inverse(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()

        @test inverse(B) === B
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
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        expected_message = "The dual of Zero is not well-defined"

        test_dim = 5
        @test_throws ErrorException(expected_message) dual(B, test_dim)
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "Zero: convert(B)" begin
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            B = Zero{precision_type_src}()
            B_converted = convert(precision_type_converted, B)
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
