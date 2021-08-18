"""
Unit tests for ==(x, y) and !=(x, y)

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

# --- File inclusions

# Test utilities
include("test_utils.jl")

# --- Tests

# ------ B::Pseudoscalar

@testset "==(B::Pseudoscalar, C::Pseudoscalar)" begin
    # Preparations
    dim = 10

    test_value = get_random_value(1)  # add 1 to keep value away from 0

    float64_or_bigfloat = (Float64, BigFloat)

    # dim(B) == dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value))
            if precision_type1 == precision_type2
                @test B == C
            elseif precision_type1 in float64_or_bigfloat &&
                    precision_type2 in float64_or_bigfloat
                @test B == C
            else
                @test B != C
            end
        end
    end

    # dim(B) != dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim + 1, precision_type2(test_value))
            @test B != C
        end
    end

    # dim(B) == dim(C), value(B) != value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value) + 1)
            @test B != C
        end
    end
end

# ------ B::Scalar

@testset "==(B::Scalar, C::Scalar)" begin
    # Preparations

    test_value = get_random_value(2)  # add 2 to keep value away from 0 and 1
    float64_or_bigfloat = (Float64, BigFloat)

    # B::Scalar, C::Scalar
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            # ==(B, C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(test_value))
            if precision_type1 == precision_type2
                @test B == C
            elseif precision_type1 in float64_or_bigfloat &&
                   precision_type2 in float64_or_bigfloat
                @test B == C
            else
                @test B != C
            end

            # value(B) != value(C)
            B = Scalar(precision_type1(test_value))
            C = Scalar(precision_type2(2 * test_value))
            @test B != C
        end
    end
end

@testset "==(B::Scalar, C::One)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Scalar, C::Zero)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Scalar{precision_type}(test_value)
        C = Zero{precision_type}()
        @test B != C
    end
end

# ------ B::One

@testset "==(B::One, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::One, C::One)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = One{precision_type1}()
            C = One{precision_type2}()
            if precision_type1 == precision_type2
                @test B === C
            else
                @test B == C
            end
        end
    end
end

@testset "==(B::One, C::Zero)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()
        C = Zero{precision_type}()
        @test B != C
    end
end

# ------ B::Zero

@testset "==(B::Zero, C::Scalar)" begin
    test_value = 5
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = Scalar{precision_type}(test_value)
        @test B != C
    end
end

@testset "==(B::Zero, C::One)" begin
    for precision_type in subtypes(AbstractFloat)
        B = Zero{precision_type}()
        C = One{precision_type}()
        @test B != C
    end
end

@testset "==(B::Zero, C::Zero)" begin
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Zero{precision_type1}()
            C = Zero{precision_type2}()
            if precision_type1 == precision_type2
                @test B === C
            else
                @test B == C
            end
        end
    end
end
