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
    test_dim = 10

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

@testset "Pseudoscalar: reverse(B)" begin
    # --- Preparations

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    for precision_type in subtypes(AbstractFloat)
        for test_dim in 5:8
            B = Pseudoscalar(test_dim, test_value)
            reverse_B = reverse(B)
            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar(test_dim, test_value) :
                Pseudoscalar(test_dim, -test_value)
            @test reverse(B) === expected_result

            @test B * reverse(B) ≈ norm(B)^2
        end
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
            expected_result = Scalar{precision_type}(precision_type(test_value))
            @test dual_B isa Scalar{precision_type}
            @test dual_B == expected_result

            # Check dual(dual_B) = (-1)^(dim_B * (dim_B - 1) / 2) B
            if mod(dim(B), 4) < 2
                @test dual(dual_B, dim=dim(B)) == B
            else
                @test dual(dual_B, dim=dim(B)) == -B
            end
        end
    end
end

@testset "Pseudoscalar: reciprocal(B)" begin
    # Preparations
    test_value = rand() + 2  # add 2 to avoid 0 and 1
    test_value = (rand() > 0.5) ? test_value : -test_value

    # Tests
    for precision_type in subtypes(AbstractFloat)
        converted_value = precision_type(test_value)

        for test_dim = 5:8
            B = Pseudoscalar(test_dim, converted_value)

            reciprocal_B = reciprocal(B)
            @test reciprocal_B isa Pseudoscalar{precision_type}

            expected_result = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1 / converted_value) :
                Pseudoscalar{precision_type}(test_dim, -1 / converted_value)
            @test reciprocal_B == expected_result

            @test B * reciprocal_B ≈ 1
        end
    end
end

@testset "Pseudoscalar: +(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

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
    test_dim = 10

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

    # Test dimension
    test_dim = 10

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

    # B::Pseudoscalar, C::Scalar
    # B::Scalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B * C == expected_result
    @test C * B == expected_result

    # B::Pseudoscalar, C::Real
    # C::Real, B::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test B * C == expected_result
    @test C * B == expected_result
end

@testset "Pseudoscalar: /(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

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
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) == dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = zero(B)
    @test wedge(B, C) == expected_result
    @test B ∧ C == expected_result

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch B ∧ C

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = Scalar(test_value_2)

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test wedge(B, C) == expected_result
    @test wedge(C, B) == expected_result
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result

    # --- B::Pseudoscalar, C::Real
    #     B::Real, C::Pseudoscalar
    B = Pseudoscalar(test_dim, test_value_1)
    C = test_value_2

    expected_result = Pseudoscalar(test_dim, test_value_1 * test_value_2)
    @test wedge(B, C) == expected_result
    @test wedge(C, B) == expected_result
    @test B ∧ C == expected_result
    @test C ∧ B == expected_result
end

@testset "Pseudoscalar: contractl(B, C), B < C" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # B::Pseudoscalar, C::Pseudoscalar

    for test_dim = 5:8
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_contractl_C = contractl(B, C)

        expected_result = mod(test_dim, 4) < 2 ?
            test_value_1 * test_value_2 :
           -test_value_1 * test_value_2

        @test B_contractl_C isa Scalar
        @test B_contractl_C == expected_result
        @test (B < C) === B_contractl_C
    end
end

@testset "Pseudoscalar: proj(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- Tests

    # --- B::Pseudoscalar, C::Pseudoscalar

    for test_dim in 5:8
        # dim(B) == dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)

        B_proj_C = proj(B, C)
        expected_result = B
        @test B_proj_C === expected_result

        # dim(B) != dim(C)
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim + 1, test_value_2)
        @test_throws DimensionMismatch proj(B, C)
    end

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    B = Scalar(test_value_1)
    C = Pseudoscalar(test_dim, test_value_2)

    expected_result = test_value_1
    B_proj_C = proj(B, C)
    @test B_proj_C isa AbstractScalar
    @test B_proj_C == expected_result

    expected_result = zero(C)
    C_proj_B = proj(C, B)
    @test C_proj_B === expected_result
end

@testset "Pseudoscalar: dual(B, C)" begin
    # --- Preparations

    # Test dimension
    test_dim = 10

    # Test values
    test_value_1 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_1 = rand() > 0.5 ? test_value_1 : -test_value_1

    test_value_2 = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value_2 = rand() > 0.5 ? test_value_2 : -test_value_2

    # --- B::Pseudoscalar, C::Pseudoscalar

    # dim(B) === dim(C)
    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(test_dim, test_value_1)
        C = Pseudoscalar(test_dim, test_value_2)
        dual_B = dual(B, C)
        expected_result = test_value_1
        @test dual_B isa Scalar
        @test dual_B == expected_result

        # Check dual(dual_B, C) = (-1)^(grade(C) * (grade(C) - 1) / 2) B
        if mod(grade(C), 4) < 2
            @test dual(dual_B, C) == B
        else
            @test dual(dual_B, C) == -B
        end
    end

    # dim(B) != dim(C)
    B = Pseudoscalar(test_dim, test_value_1)
    C = Pseudoscalar(test_dim + 1, test_value_2)
    @test_throws DimensionMismatch dual(B, C)

    # --- B::Pseudoscalar, C::Scalar
    #     B::Scalar, C::Pseudoscalar

    C = Scalar(test_value_2)

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result
    end

    # --- B::Pseudoscalar, C::One
    #     B::One, C::Pseudoscalar

    C = One()

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result
    end

    # --- B::Pseudoscalar, C::Real
    #     B::Real, C::Pseudoscalar

    C = test_value_2

    for dim_B in test_dim:test_dim + 3
        B = Pseudoscalar(dim_B, test_value_1)
        expected_result = zero(B)
        @test dual(B, C) == expected_result
    end
end
