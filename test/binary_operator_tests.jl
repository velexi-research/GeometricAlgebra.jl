"""
Unit tests for AbstractBlade binary operators.

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
using Test

# GeometricAlgebra.jl
using GeometricAlgebra


# --- ∧(B, C)

@testset "∧(B, C) tests" begin
    # --- B, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vectors = hcat([0; 0; 3; 0; 0],
                     [0; 0; 0; 4; 0])
    C = Blade(C_vectors)

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C ≈ Blade(hcat(B_vectors, C_vectors))

    # --- B::Blade, C::Scalar
    #     B::Scalar, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value
    C = Scalar(value)

    expected_B_wedge_C = Blade(B_vectors, volume=2 * value)

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C ≈ expected_B_wedge_C

    C_wedge_B = C ∧ B
    @test C_wedge_B == outer(C, B)
    @test C_wedge_B ≈ expected_B_wedge_C

    # --- B::Blade, C::Zero
    #     B::Zero, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C = Zero()

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C === Zero()

    C_wedge_B = C ∧ B
    @test C_wedge_B == outer(C, B)
    @test C_wedge_B === Zero()

    # --- B::Blade, C::Vector
    #     B::Vector, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C_vector = [0; 0; 3; 0; 0]
    C = Blade(C_vector)

    expected_B_wedge_C = Blade(hcat(B_vectors, C_vector))

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C ≈ expected_B_wedge_C

    C_wedge_B = C ∧ B
    @test C_wedge_B == outer(C, B)
    @test C_wedge_B ≈ (-1)^(grade(B)) * expected_B_wedge_C

    # --- B::Blade, C::Real
    #     B::Real, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    C = 2.5

    expected_B_wedge_C = Blade(B_vectors, volume=5.0)

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C == expected_B_wedge_C

    C_wedge_B = C ∧ B
    @test C_wedge_B == outer(C, B)
    @test C_wedge_B == expected_B_wedge_C
end

# --- *(B, C)

@testset "*(x, B) tests: x::Real, B::Blade" begin
    # Preparations
    x = rand() + 1  # add 1 to avoid 0
    x = rand() > 0.5 ? x : -x

    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    expected_x_times_B = Blade(B_vectors, volume=2 * x)

    # Exercise functionality
    @test x * B == expected_x_times_B
    @test B * x == expected_x_times_B
end

@testset "*(S, B) tests: x::Scalar, B::Blade" begin
    # --- Preparations

    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    # --- S::Scalar, B::Blade
    #     B::Blade, S::Scalar

    # Preparations
    x = rand() + 1  # add 1 to avoid 0
    x = rand() > 0.5 ? x : -x
    S = Scalar(x)

    expected_S_times_B = Blade(B_vectors, volume=2 * x)

    # Exercise functionality
    @test S * B == expected_S_times_B
    @test B * S == expected_S_times_B

    # --- S::Zero, B::Blade
    #     B::Blade, S::Zero

    @test Zero() * B === Zero()
    @test B * Zero() === Zero()
end
