#   Copyright 2020 Velexi Corporation
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
Unit tests for the One type.
"""

# --- Imports

# Standard library
import InteractiveUtils.subtypes
using Test

# GeometricAlgebra.jl
using GeometricAlgebra

# --- Constructor tests

@testset "One: constructor tests" begin
    # --- Inner constructor

    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        @test B isa One{precision_type}
        @test B === One{precision_type}()
    end

    # --- Outer constructor

    B = One()
    @test B isa One{Float64}
end

@testset "One: one()" begin

    # --- one(M::AbstractMultivector)

    for precision_type in subtypes(AbstractFloat)
        # --- M::AbstractScalar

        # M::Zero
        B = one(Zero{precision_type}())
        @test B isa One{precision_type}

        # M::One
        B = one(One{precision_type}())
        @test B isa One{precision_type}

        # M::Scalar
        B = one(Scalar{precision_type}(5))
        @test B isa One{precision_type}

        # --- M::AbstractBlade

        # M::Blade
        B = one(Blade{precision_type}([1 2 3]))
        @test B isa One{precision_type}

        # M::Pseudoscalar
        B = one(Pseudoscalar{precision_type}(10, 5))
        @test B isa One{precision_type}

        # --- M::AbstractMultivector

        # M::Multivector
        B = one(Multivector{precision_type}([Scalar(3)]))
        @test B isa One{precision_type}
    end

    # --- one(M::Type{<:AbstractMultivector})

    B = one(AbstractScalar)
    @test B isa One{Float64}

    B = one(Zero)
    @test B isa One{Float64}

    B = one(One)
    @test B isa One{Float64}

    B = one(Scalar)
    @test B isa One{Float64}

    B = one(AbstractBlade)
    @test B isa One{Float64}

    B = one(Blade)
    @test B isa One{Float64}

    B = one(Pseudoscalar)
    @test B isa One{Float64}

    B = one(AbstractMultivector)
    @test B isa One{Float64}

    B = one(Multivector)
    @test B isa One{Float64}

    # --- one(::Type{<:AbstractMultivector{T}})

    for precision_type in subtypes(AbstractFloat)
        # --- Type{<:AbstractScalar{T}}

        B = one(AbstractScalar{precision_type})
        @test B isa One{precision_type}

        B = one(Zero{precision_type})
        @test B isa One{precision_type}

        B = one(One{precision_type})
        @test B isa One{precision_type}

        B = one(Scalar{precision_type})
        @test B isa One{precision_type}

        # --- Type{<:AbstractBlade{T}}

        B = one(AbstractBlade{precision_type})
        @test B isa One{precision_type}

        B = one(Blade{precision_type})
        @test B isa One{precision_type}

        B = one(Pseudoscalar{precision_type})
        @test B isa One{precision_type}

        # --- Type{<:AbstractMultivector{T}}

        B = one(AbstractMultivector{precision_type})
        @test B isa One{precision_type}

        B = one(Multivector{precision_type})
        @test B isa One{precision_type}
    end
end

# --- Test comparison operations

@testset "One: isone(B)" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        # Blade
        vectors = [3 3; 4 4; 0 1; 0 1]
        B = Blade{precision_type}(vectors)
        @test !isone(B)

        # Pseudoscalar
        B = Pseudoscalar{precision_type}(10, 3)
        @test !isone(B)

        # Scalar
        B = Scalar{precision_type}(1)
        @test isone(B)

        B = Scalar{precision_type}(3)
        @test !isone(B)

        # Zero
        B = Zero{precision_type}()
        @test !isone(B)

        # One
        B = One{precision_type}()
        @test isone(B)
    end
end

# --- Test attribute methods

@testset "One: AbstractMultivector attribute functions" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        @test dim(B) == 0
        @test grades(B) == [0]
        @test blades(B) == [B]
        @test norm(B) isa precision_type
        @test norm(B) == 1

        @test B[0] == [B]
        @test B[1] == []
    end
end

@testset "One: AbstractBlade attribute functions" begin
    # Basic functions
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        @test grade(B) == 0
        @test basis(B) == 1
        @test volume(B) isa precision_type
        @test volume(B) == 1
        @test sign(B) == 1
    end
end

@testset "One: AbstractScalar attribute functions" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        @test value(B) isa precision_type
        @test value(B) == 1
    end
end

# --- Tests for AbstractMultivector interface functions

@testset "One: -(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        negative_B = -B
        @test negative_B isa Scalar{precision_type}
        @test negative_B == Scalar{precision_type}(-1)
    end
end

@testset "One: reverse(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        reverse_B = reverse(B)
        @test reverse_B === B
        @test B * reverse_B == norm(B)^2
    end
end

@testset "One: dual(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        for test_dim in 5:8
            dual_B = dual(B, test_dim)

            expected_dual = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1) :
                Pseudoscalar{precision_type}(test_dim, -1)
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_dual
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "One: convert(B)" begin
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            B = One{precision_type_src}()
            B_converted = convert(precision_type_converted, B)
            @test B_converted isa One{precision_type_converted}
        end
    end
end

# --- Tests for AbstractBlade interface functions

@testset "One: inv(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        inverse_B = inv(B)
        @test inverse_B === B
    end
end
