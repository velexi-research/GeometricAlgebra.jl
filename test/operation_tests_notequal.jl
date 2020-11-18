"""
Unit tests for !=(x, y)

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

# ------ M::Multivector

@testset "!=(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "!=(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "!=(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Blade, C::Pseudoscalar)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade{precision_type1}(vectors)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
        end
    end
end

@testset "!=(B::Blade, C::Scalar)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Blade{precision_type1}(vectors)
            C = Scalar{precision_type2}(test_value)
            @test B != C
        end
    end
end

@testset "!=(B::Blade, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Blade, C::Zero)" begin
    @test_skip 1
end

@testset "!=(B::Blade, C::Real)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(vectors)
        C = test_value
        @test B != C
    end
end

@testset "!=(B::Blade, C::Vector)" begin
    @test_skip 1
end

# ------ B::Pseudoscalar

@testset "!=(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Pseudoscalar, C::Blade)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar{precision_type2}(test_dim, test_value)
            C = Blade{precision_type1}(vectors)
            @test B != C
        end
    end
end

@testset "!=(B::Pseudoscalar, C::Scalar)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar{precision_type2}(test_dim, test_value)
            C = Scalar{precision_type1}(test_value)
            @test B != C
        end
    end
end

@testset "!=(B::Pseudoscalar, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Pseudoscalar, C::Zero)" begin
    @test_skip 1
end

@testset "!=(B::Pseudoscalar, C::Real)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = test_value
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::Scalar

@testset "!=(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Scalar, C::Blade)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Scalar{precision_type2}(test_value)
            C = Blade{precision_type1}(vectors)
            @test B != C
        end
    end
end

@testset "!=(B::Scalar, C::Pseudoscalar)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Scalar{precision_type1}(test_value)
            C = Pseudoscalar{precision_type2}(test_dim, test_value)
            @test B != C
        end
    end
end

@testset "!=(B::Scalar, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Scalar, C::Zero)" begin
    @test_skip 1
end

@testset "!=(B::Scalar, C::Real)" begin
    @test_skip 1
end

@testset "!=(B::Scalar, C::Vector)" begin
    @test_skip 1
end

# ------ B::One

@testset "!=(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Blade)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Scalar)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Zero)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Real)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Vector)" begin
    @test_skip 1
end

# ------ B::Zero

@testset "!=(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Blade)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Scalar)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Real)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Vector)" begin
    @test_skip 1
end

# ------ B::Real

@testset "!=(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Real, C::Blade)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = test_value
        C = Blade{precision_type}(vectors)
        @test B != C
    end
end

@testset "!=(B::Real, C::Pseudoscalar)" begin
    vectors = [3 3; 4 4; 0 1]
    test_dim = size(vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = test_value
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Real, C::Scalar)" begin
    @test_skip 1
end

@testset "!=(B::Real, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Real, C::Zero)" begin
    @test_skip 1
end

# ------ B::Vector

@testset "!=(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::Blade)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::Scalar)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::One)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::Zero)" begin
    @test_skip 1
end
