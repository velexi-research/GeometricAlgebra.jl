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

    # --- B::Blade, C::AbstractScalar
    #     B::AbstractScalar, C::Blade

    # Preparations
    B_vectors = hcat([1; 0; 0; 0; 0],
                     [0; 2; 0; 0; 0])
    B = Blade(B_vectors)

    value = rand() + 1  # add 1 to avoid 0
    value = rand() > 0.5 ? value : -value
    C = Scalar(value)

    # Exercise functionality
    B_wedge_C = B ∧ C
    @test B_wedge_C == outer(B, C)
    @test B_wedge_C ≈ Blade(B_vectors, volume=2 * value)
end
