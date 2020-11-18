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
    test_vectors = [3 3; 4 4; 0 1]

    test_dim = size(test_vectors, 1)
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::Scalar)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::One)" begin
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Blade, C::Zero)" begin
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Blade, C::Real)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "!=(B::Blade, C::Vector)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_vector = [1; 2; 3]

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Pseudoscalar

@testset "!=(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Pseudoscalar, C::Blade)" begin
    test_dim = 3
    test_value = 5

    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::Scalar)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::One)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::Zero)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::Real)" begin
    test_dim = 10
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "!=(B::Pseudoscalar, C::Vector)" begin
    test_dim = 3
    test_value = 5

    test_vector = [3; 4; 1]

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Scalar

@testset "!=(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Scalar, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::Scalar, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10

    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Scalar, C::One)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Scalar, C::Zero)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Scalar, C::Real)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = precision_type(test_value + 1)
        @test B != C
    end
end

@testset "!=(B::Scalar, C::Vector)" begin
    test_value = 5
    test_vector = [3; 4; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::One

@testset "!=(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::One, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::One, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::One, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::One, C::Zero)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Zero{precision_type}()
        @test B != C
    end
end

@testset "!=(B::One, C::Real)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "!=(B::One, C::Vector)" begin
    test_vector = [3; 4; 0]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Zero

@testset "!=(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Zero, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::Zero, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Zero, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::Zero, C::One)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Zero, C::Real)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = precision_type(test_value)
        @test B != C
    end
end

@testset "!=(B::Zero, C::Vector)" begin
    test_vector = [3; 4; 0]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Vector{precision_type}(test_vector)
        @test B != C
    end
end

# ------ B::Real

@testset "!=(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Real, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::Real, C::Pseudoscalar)" begin
    test_value = 5
    test_dim = 10
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Real, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Scalar{precision_type}(test_value + 1)
        @test B != C
    end
end

@testset "!=(B::Real, C::One)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Real, C::Zero)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

# ------ B::Vector

@testset "!=(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "!=(B::Vector, C::Blade)" begin
    test_vector = [1; 2; 3]
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Blade{precision_type}(test_vectors)
        @test B != C
    end
end

@testset "!=(B::Vector, C::Pseudoscalar)" begin
    test_vector = [1; 2; 3]

    test_dim = 3
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B != C
    end
end

@testset "!=(B::Vector, C::Scalar)" begin
    test_vector = [1; 2; 3]
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "!=(B::Vector, C::One)" begin
    test_vector = [1; 2; 3]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "!=(B::Vector, C::Zero)" begin
    test_vector = [1; 2; 3]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Zero{precision_type}()
        @test B != C
    end
end
