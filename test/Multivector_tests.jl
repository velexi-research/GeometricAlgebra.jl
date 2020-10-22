"""
Unit tests for the Multivector type.

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

# External packages
import DataStructures.SDMKeyIteration, DataStructures.SortedDict

# GeometricAlgebra.jl
using GeometricAlgebra


# --- Constructor tests

@testset "Multivector: inner constructor" begin
    # Notes
    # -----
    # * Test value of constructed instance

    # --- Preparations

    vectors = [3 3; -4 -4; 0 1]
    one_vector = [3; 4; 0]

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = length(one_vector)

    # --- Multivector{T}(blades::Vector{AbstractBlade{T}};
    #                    reduced::Bool=false) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        scalar = Scalar{precision_type}(test_value)
        one_blade = Blade{precision_type}(one_vector)
        two_blade = Blade{precision_type}(vectors)
        pseudoscalar = Pseudoscalar{precision_type}(test_dim, test_value)

        blades = Vector([scalar, one_blade, two_blade, pseudoscalar])

        # default value for `reduced`
        M = Multivector{precision_type}(blades)
        @test length(M.parts) == length(blades)
        for B in blades
            @test haskey(M.parts, grade(B))
            @test length(M.parts[grade(B)]) == 1
        end
    end
end

@testset "Multivector: outer constructor" begin
    # Notes
    # -----
    # * Test type of constructed instances. Correct construction of instances
    #   is tested by the inner constructor tests.

    # --- Preparations

    vectors = [3 3; -4 -4; 0 1]
    one_vector = [3; 4; 0]

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = length(one_vector)

    # --- Multivector{T}(blades::Vector{AbstractBlade{T}};
    #                    reduced::Bool=false) where {T<:AbstractFloat}

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        scalar = Scalar{precision_type}(test_value)
        two_blade = Blade{precision_type}(vectors)
        pseudoscalar = Pseudoscalar{precision_type}(test_dim, test_value)

        blades = Vector([scalar, two_blade, pseudoscalar])

        # default value for `reduced`
        M = Multivector(blades)
        @test M isa Multivector{precision_type}
    end
end

# --- Function tests

@testset "AbstractMultivector interface: M::Multivector" begin
    # --- Preparations

    vectors = [3 3; -4 -4; 0 1]
    one_vector = [3; 4; 0]

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = length(one_vector)

    # --- Test basic functions

    for precision_type in subtypes(AbstractFloat)
        # Preparations
        scalar = Scalar{precision_type}(test_value)
        one_blade = Blade{precision_type}(one_vector)
        pseudoscalar = Pseudoscalar{precision_type}(test_dim, test_value)
        expected_blades = Vector([scalar, one_blade, pseudoscalar])
        M = Multivector{precision_type}(expected_blades)

        expected_parts = SortedDict(0=>Vector([scalar]),
                                    1=>Vector([one_blade]),
                                    test_dim=>Vector([pseudoscalar]))
        expected_blades = reduce(vcat, values(expected_parts))

        # grades()
        @test grades(M) == [0, 1, 3]
        @test grades(M, collect=false) isa SDMKeyIteration

        # blades()
        @test blades(M) isa Vector{AbstractBlade}
        @test blades(M) == expected_blades

        # getindex()
        for k in 1:test_dim
            if k in grades(M)
                k_vectors = M[k]
                @test k_vectors isa Vector{AbstractBlade}
                @test k_vectors == expected_parts[k]
            end
        end

        # norm()
        @test_skip norm(M) == 0
    end
end

@testset "convert(M): M::Multivector" begin
    # Preparations
    vectors = [3 3; -4 -4; 0 1]
    one_vector = [3; 4; 0]

    test_value = rand() + 1  # add 1 to avoid 0
    test_value = rand() > 0.5 ? test_value : -test_value

    test_dim = length(one_vector)

    # Tests
    for precision_type_converted in subtypes(AbstractFloat)
        for precision_type_src in subtypes(AbstractFloat)
            # Preparations
            scalar = Scalar{precision_type_src}(test_value)
            one_blade = Blade{precision_type_src}(one_vector)
            pseudoscalar = Pseudoscalar{precision_type_src}(test_dim,
                                                            test_value)
            blades = Vector([scalar, one_blade, pseudoscalar])
            M = Multivector{precision_type_src}(blades)

            # Exercise functionality and check results
            M_converted = convert(Multivector{precision_type_converted}, M)

            @test M_converted isa Multivector{precision_type_converted}

            if precision_type_src == precision_type_converted
                @test M_converted === M
            else
                @test M_converted !== M
                @test_skip M_converted ≈ M
            end
        end
    end
end
