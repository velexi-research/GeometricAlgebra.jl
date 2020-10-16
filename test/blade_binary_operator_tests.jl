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


# --- *(B, C)

@testset "*(x, B): x::Real, B::Blade" begin
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

@testset "*(B, C): B::Blade, C::Scalar" begin
    # --- Preparations

    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    x = rand() + 1  # add 1 to avoid 0
    x = rand() > 0.5 ? x : -x
    C = Scalar(x)

    expected_B_times_C = Blade(B_vectors, volume=2 * x)

    # Exercise functionality
    @test B * C == expected_B_times_C
    @test C * B == expected_B_times_C
end
