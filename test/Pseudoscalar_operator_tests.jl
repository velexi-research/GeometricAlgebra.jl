"""
Operator unit tests for Pseudoscalar type.

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

@testset "Pseudoscalar: ==(B, C)" begin
    # Preparations
    dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

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

@testset "Pseudoscalar: ≈(B, C)" begin
    # Preparations
    dim = 10

    test_value = rand()
    test_value = rand() > 0.5 ? test_value : -test_value

    # dim(B) == dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value))
            @test B ≈ C
        end
    end

    # dim(B) != dim(C), value(B) == value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim + 1, precision_type2(test_value))
            @test B ≉ C
        end
    end

    # dim(B) == dim(C), value(B) != value(C)
    for precision_type1 in subtypes(AbstractFloat)
        for precision_type2 in subtypes(AbstractFloat)
            B = Pseudoscalar(dim, precision_type1(test_value))
            C = Pseudoscalar(dim, precision_type2(test_value) + 1)
            @test B ≉ C
        end
    end
end

@testset "Pseudoscalar: -(B)" begin
    # Preparations
    test_dim = 5

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)

        minus_B = -B
        @test minus_B isa Pseudoscalar{precision_type}
        @test minus_B == Pseudoscalar{precision_type}(test_dim, -test_value)
    end
end

@testset "Pseudoscalar: dual(B)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        for test_dim in 5:8
            B = Pseudoscalar{precision_type}(test_dim, test_value)

            dual_B = dual(B)
            expected_dual = Scalar{precision_type}(precision_type(test_value))
            @test dual_B isa Scalar{precision_type}
            @test dual_B == expected_dual
        end
    end
end

@testset "Pseudoscalar: reciprocal(B)" begin
    # Preparations
    test_dim = 5

    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        B = Pseudoscalar{precision_type}(test_dim, test_value)

        reciprocal_B = reciprocal(B)
        @test reciprocal_B isa Pseudoscalar{precision_type}
        @test reciprocal_B ==
            Pseudoscalar{precision_type}(test_dim,
                                         1 / precision_type(test_value))
    end
end

@testset "Pseudoscalar: +(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 5

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_plus_C = B + C
    expected_result = Pseudoscalar(test_dim, test_value_1 + test_value_2)
    @test B_plus_C isa Pseudoscalar
    @test B_plus_C == expected_result
end

@testset "Pseudoscalar: -(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 5

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = abs(test_value_1) + rand() + 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # ------ B::Pseudoscalar, C::Pseudoscalar

    # B != C
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_minus_C = B - C
    expected_result = Pseudoscalar(test_dim, test_value_1 - test_value_2)
    @test B_minus_C isa Pseudoscalar
    @test B_minus_C == expected_result

    # B == C
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_1)

    B_minus_C = B - C
    expected_result = Zero()
    @test B_minus_C === expected_result
end

@testset "Pseudoscalar: *(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_times_C = B * C
        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_times_C isa Scalar
        @test B_times_C == expected_result
    end
end

@testset "Pseudoscalar: /(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 5

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_slash_C = B / C
    expected_result = Scalar(test_value_1 / test_value_2)
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result
end

@testset "Pseudoscalar: wedge(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 5

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_wedge_C = wedge(B, C)
    expected_result = zero(B)
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C
end

@testset "Pseudoscalar: contractl(B, C), B < C" begin
    # --- Preparations

    # Test dimension
    test_dim = 5

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    B_contractl_C = contractl(B, C)
    expected_result = test_value_1 * test_value_2
    @test B_contractl_C isa Scalar
    @test B_contractl_C == expected_result
    @test (B < C) === B_contractl_C
end

@testset "Pseudoscalar: proj(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_proj_C = proj(B, C)
        expected_result = B
        @test B_proj_C === expected_result
    end
end

@testset "Pseudoscalar: dual(B, C)" begin
    # --- Preparations

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # C::Pseudoscalar
    for test_dim in 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)
        B_dual_C = dual(B, C)
        expected_result = test_value_1
        @test B_dual_C isa Scalar
        @test B_dual_C == expected_result
    end
end
