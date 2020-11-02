"""
Unit tests for !=(M, N)

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

@testset "!=(M, N)" begin
    # Preparations
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- B::AbstractBlade, C::AbstractBlade

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # B::Blade, C::Scalar
            # B::Scalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Scalar{precision_type2}(test_value)
            @test B != C
            @test C != B

            # B::Blade, C::Pseudoscalar
            # B::Pseudoscalar, C::Blade
            B = Blade{precision_type1}(vectors)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
            @test C != B

            # B::Scalar, C::Pseudoscalar
            # B::Pseudoscalar, C::Scalar
            B = Scalar{precision_type1}(test_value)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
            @test C != B
        end
    end

    # --- B::Union{Blade, Pseudoscalar}, C::Real
    #     B::Real, C::Union{Blade, Pseudoscalar}

    for precision_type in subtypes(AbstractFloat)
        # B::Blade, C::Real
        # B::Real, C::Blade
        B = Blade{precision_type}(vectors)
        @test B != test_value
        @test test_value != B

        # B::Pseudoscalar, C::Real
        # B::Real, C::Pseudoscalar
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != test_value
        @test test_value != B
    end
end
