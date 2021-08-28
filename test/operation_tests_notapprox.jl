"""
Unit tests for !isapprox(x, y)

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

#=
# ------ M::Multivector

@testset "!isapprox(M::Multivector, N::Blade)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Pseudoscalar)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Scalar)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::One)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Zero)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Real)" begin
    @test_skip 1
end

@testset "!isapprox(M::Multivector, N::Vector)" begin
    @test_skip 1
end

# ------ B::Blade

@testset "!isapprox(B::Blade, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Blade, C::Pseudoscalar)" begin
    test_vectors = [3 3; 4 4; 0 1]

    test_dim = 3
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Pseudoscalar{precision_type}(test_dim, test_value)
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Scalar)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5

    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Scalar{precision_type}(test_value)
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::One)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = One{precision_type}()
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Zero)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Zero{precision_type}()
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Real)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = precision_type(test_value)
        @test B ≉ C
    end
end

@testset "!isapprox(B::Blade, C::Vector)" begin
    test_vectors = [3 3; 4 4; 0 1]
    test_vector = [3; 4; 1]
    for precision_type in subtypes(AbstractFloat)
        B = Blade{precision_type}(test_vectors)
        C = Vector{precision_type}(test_vector)
        @test B ≉ C
    end
end

# ------ B::Pseudoscalar

@testset "!isapprox(B::Pseudoscalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Pseudoscalar, C::Blade)" begin
    test_dim = 3
    test_value = 5

    test_vectors = [3 3; 4 4; 0 1]

    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

# ------ B::Scalar

@testset "!isapprox(B::Scalar, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Scalar, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

# ------ B::One

@testset "!isapprox(B::One, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::One, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

# ------ B::Zero

@testset "!isapprox(B::Zero, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Zero, C::Blade)" begin
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

# ------ B::Real

@testset "!isapprox(B::Real, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Real, C::Blade)" begin
    test_value = 5
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = precision_type(test_value)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end

# ------ B::Vector

@testset "!isapprox(B::Vector, M::Multivector)" begin
    @test_skip 1
end

@testset "!isapprox(B::Vector, C::Blade)" begin
    test_vector = [3; 4; 1]
    test_vectors = [3 3; 4 4; 0 1]
    for precision_type in subtypes(AbstractFloat)
        B = Vector{precision_type}(test_vector)
        C = Blade{precision_type}(test_vectors)
        @test B ≉ C
    end
end
=#
