"""
Operator unit tests for One type.

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

@testset "-(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        minus_B = -B
        @test minus_B isa Scalar{precision_type}
        @test minus_B == Scalar{precision_type}(-1)
    end
end

@testset "dual(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        for test_dim in 5:8
            dual_B = dual(B, dim=test_dim)

            expected_dual = mod(test_dim, 4) < 2 ?
                Pseudoscalar{precision_type}(test_dim, 1) :
                Pseudoscalar{precision_type}(test_dim, -1)
            @test dual_B isa Pseudoscalar{precision_type}
            @test dual_B == expected_dual
        end

        expected_message = "The dual of a scalar is not well-defined if " *
                           "`dim` is not specified"
        @test_throws ErrorException(expected_message) dual(B)
    end
end

@testset "reciprocal(B)" begin
    for precision_type in subtypes(AbstractFloat)
        B = One{precision_type}()

        reciprocal_B = reciprocal(B)
        @test reciprocal_B === B
    end
end

@testset "+(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    expected_result = 2
    C_plus_B = C + B
    @test B + C isa Scalar
    @test B + C == expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_plus_C = B + C
    expected_result = 1 + test_value
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_plus_C = B + C
    expected_result = 1 + test_value
    @test B_plus_C isa Scalar
    @test B_plus_C == expected_result

    C_plus_B = C + B
    expected_result = 1 + test_value
    @test C_plus_B isa Scalar
    @test C_plus_B == expected_result
end

@testset "-(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # ------ B::One, C::One

    B = One()
    C = One()

    expected_result = Zero()
    @test B - C === expected_result

    # ------ B::One, C::Scalar
    #        B::Scalar, C::One

    B = One()
    C = Scalar(test_value)

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # ------ B::One, C::Real
    #        B::Real, C::One

    # B != C
    B = One()
    C = test_value

    B_minus_C = B - C
    expected_result = 1 - test_value
    @test B_minus_C isa Scalar
    @test B_minus_C == expected_result

    C_minus_B = C - B
    expected_result = test_value - 1
    @test C_minus_B isa Scalar
    @test C_minus_B == expected_result

    # B == C
    B = One()
    C = 1

    expected_result = Zero()

    B_minus_C = B - C
    @test B_minus_C === expected_result

    C_minus_B = C - B
    @test C_minus_B === expected_result
end

@testset "*(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    expected_result = One()
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    expected_result = C
    @test B * C === expected_result
    @test C * B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_times_C = B * C
    expected_result = test_value
    @test B_times_C isa Scalar
    @test B_times_C == expected_result

    C_times_B = C * B
    expected_result = test_value
    @test C_times_B isa Scalar
    @test C_times_B == expected_result
end

@testset "/(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_slash_C = B / C
    expected_result = One()
    @test B_slash_C === expected_result

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_slash_C = B / C
    expected_result = 1 / test_value
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_slash_C = B / C
    expected_result = 1 / test_value
    @test B_slash_C isa Scalar
    @test B_slash_C == expected_result

    C_slash_B = C / B
    expected_result = test_value
    @test C_slash_B isa Scalar
    @test C_slash_B == expected_result
end

@testset "wedge(B, C), B ∧ C" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_wedge_C = wedge(B, C)
    expected_result = One()
    @test B_wedge_C === expected_result
    @test B ∧ C === B_wedge_C

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_wedge_C = wedge(B, C)
    expected_result = C
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = C
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B === C_wedge_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_wedge_C = wedge(B, C)
    expected_result = C
    @test B_wedge_C isa Scalar
    @test B_wedge_C == expected_result
    @test B ∧ C === B_wedge_C

    C_wedge_B = wedge(C, B)
    expected_result = C
    @test C_wedge_B isa Scalar
    @test C_wedge_B == expected_result
    @test C ∧ B === C_wedge_B
end

@testset "contractl(B, C), B < C" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_left_contract_C = contractl(B, C)
    expected_result = One()
    @test B_left_contract_C === expected_result
    @test (B < C) === B_left_contract_C

    # B::One, C::Scalar
    # B::Scalar, C::One
    B = One()
    C = Scalar(test_value)

    B_left_contract_C = contractl(B, C)
    expected_result = C
    @test B_left_contract_C isa Scalar
    @test B_left_contract_C == expected_result
    @test (B < C) === B_left_contract_C

    C_left_contract_B = contractl(C, B)
    expected_result = C
    @test C_left_contract_B isa Scalar
    @test C_left_contract_B == expected_result
    @test (C < B) === C_left_contract_B

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_left_contract_C = contractl(B, C)
    expected_result = C
    @test B_left_contract_C isa Scalar
    @test B_left_contract_C == expected_result
    @test (B < C) === B_left_contract_C

    C_left_contract_B = contractl(C, B)
    expected_result = C
    @test C_left_contract_B isa Scalar
    @test C_left_contract_B == expected_result
    @test (C < B) === C_left_contract_B
end

@testset "proj(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    # --- Tests

    # B::One, C::One
    B = One()
    C = One()

    B_proj_C = proj(B, C)
    expected_result = One()
    @test B_proj_C === expected_result

    # B::Scalar, C::One
    # B::One, C::Scalar
    B = One()
    C = Scalar(test_value)

    B_proj_C = proj(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = C
    @test C_proj_B === expected_result

    # B::One, C::Real
    # B::Real, C::One
    B = One()
    C = test_value

    B_proj_C = proj(B, C)
    expected_result = B
    @test B_proj_C === expected_result

    C_proj_B = proj(C, B)
    expected_result = C
    @test C_proj_B isa Scalar
    @test C_proj_B == expected_result
end

@testset "dual(B, C)" begin
    # --- Preparations

    # Test values
    test_value = rand() + 2  # add 2 to keep value away from 0 and 1
    test_value = rand() > 0.5 ? test_value : -test_value

    B = One()

    # C::One
    C = One()
    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    # C::Scalar
    C = Scalar(test_value)
    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result

    # C::Blade
    test_dim = 10
    for test_grade in 5:8
        C = Blade(randn(test_dim, test_grade))
        B_dual_C = dual(B, C)
        expected_result = mod(grade(C), 4) < 2 ?
            Blade(C, volume=1, copy_basis=false) :
            Blade(C, volume=-1, copy_basis=false)
        @test B_dual_C isa Blade
        @test B_dual_C == expected_result
    end

    # C::Pseudoscalar
    for test_dim in 5:8
        C = Pseudoscalar(test_dim, test_value)
        B_dual_C = dual(B, C)
        expected_result = mod(grade(C), 4) < 2 ?
            Pseudoscalar(test_dim, 1) :
            Pseudoscalar(test_dim, -1)
        @test B_dual_C isa Pseudoscalar
        @test B_dual_C == expected_result
    end

    # C::Real
    C = test_value
    B_dual_C = dual(B, C)
    expected_result = B
    @test B_dual_C === expected_result
end
