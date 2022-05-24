#   Copyright 2020 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
multiply.jl defines methods for the *(x, y) function
"""

# --- Exports

import Base.:(*)

# --- Method definitions

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, N::AbstractMultivector
"""
    *(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the geometric product of `M` and `N`.
"""
*(M::AbstractMultivector, N::AbstractMultivector) = nothing # TODO: implement this method

# M::AbstractMultivector, B::AbstractBlade
# B::AbstractBlade, M::AbstractMultivector
*(M::AbstractMultivector, B::AbstractBlade) =
    Multivector(map(C -> C * B, blades(M)))

*(B::AbstractBlade, M::AbstractMultivector) =
    Multivector(map(C -> B * C, blades(M)))

# M::AbstractMultivector, B::AbstractScalar,
# B::AbstractScalar, M::AbstractMultivector
*(M::AbstractMultivector, B::AbstractScalar) = B * M
*(B::AbstractScalar, M::AbstractMultivector) = value(B) * M

# M::AbstractMultivector, B::One
# B::One, M::AbstractMultivector
*(M::AbstractMultivector, B::One) = M
*(B::One, M::AbstractMultivector) = M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
*(M::AbstractMultivector, B::Zero) = B
*(B::Zero, M::AbstractMultivector) = B

# M::AbstractMultivector, x::Real
# x::Real, M::AbstractMultivector
*(M::AbstractMultivector, x::Real) = x * M
*(x::Real, M::AbstractMultivector) = Multivector(map(B -> x * B, blades(M)))

# M::AbstractMultivector, v::Vector
# v::Vector, M::AbstractMultivector
*(M::AbstractMultivector, v::Vector) = Multivector(map(B -> B * v, blades(M)))
*(v::Vector, M::AbstractMultivector) = Multivector(map(B -> v * B, blades(M)))

# ------ Specializations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, B::AbstractBlade
*(B::AbstractBlade, C::AbstractScalar) =
    Blade(B, volume=volume(B) * value(C), copy_basis=false)
*(C::AbstractScalar, B::AbstractBlade) = B * C

# B::AbstractBlade, C::One
# B::One, C::AbstractBlade
*(B::AbstractBlade, C::One) = B
*(B::One, C::AbstractBlade) = C

# B::AbstractBlade, C::Zero
# B::Zero, C::AbstractBlade
*(B::AbstractBlade, C::Zero) = C
*(B::Zero, C::AbstractBlade) = B

# B::AbstractBlade, x::Real
# x::Real, B::AbstractBlade
*(B::AbstractBlade, x::Real) = Blade(B, volume=volume(B) * x)
*(x::Real, B::AbstractBlade) = B * x

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function *(B::Blade, C::Blade)
    # Check arguments
    assert_dim_equal(B, C)

    # Compute geometric product
    M = C
    for i in grade(B):-1:1
        M = basis(B)[:, i] * M
    end
    M = volume(B) * M
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
*(B::Blade, C::Pseudoscalar) = B ⋅ C
*(B::Pseudoscalar, C::Blade) = zero(B)

# B::Blade, v::Vector
# v::Vector, B::Blade
function *(B::Blade, v::Vector{<:Real})
    # --- Check arguments

    # Check that B and v are from the same real vector space
    assert_dim_equal(v, B)

    # Check that grade(B) == 1. B * v is not defined when grade(B) > 1
    if grade(B) > 1
        error("Geometric product B * v is not defined when grade(B) > 1")
    end

    # Compute geometric product
    B_dot_v = B ⋅ v
    B_wedge_v = B ∧ v

    if B_dot_v == zero(B_dot_v)
        return B_wedge_v
    elseif B_wedge_v == zero(B_wedge_v)
        return B_dot_v
    end

    Multivector([B_dot_v, B_wedge_v])
end

function *(v::Vector{<:Real}, B::Blade)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute geometric product
    v_dot_B = v ⋅ B
    v_wedge_B = v ∧ B

    if v_dot_B == zero(v_dot_B)
        return v_wedge_B
    elseif v_wedge_B == zero(v_wedge_B)
        return v_dot_B
    end

    Multivector([v_dot_B, v_wedge_B])
end

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
*(B::Pseudoscalar, C::Pseudoscalar) = contract_left(B, C)

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, C::Pseudoscalar
*(B::Pseudoscalar, C::AbstractScalar) =
    Pseudoscalar(B, value=value(B) * value(C))

*(B::AbstractScalar, C::Pseudoscalar) = C * B

# B::Pseudoscalar, C::One
# B::One, C::Pseudoscalar
*(B::Pseudoscalar, C::One) = B
*(B::One, C::Pseudoscalar) = C

# B::Pseudoscalar, C::Zero
# B::Zero, C::Pseudoscalar
*(B::Pseudoscalar, C::Zero) = C
*(B::Zero, C::Pseudoscalar) = B

# B::Pseudoscalar, x::Real
# x::Real, B::Pseudoscalar
*(B::Pseudoscalar, x::Real) = Pseudoscalar(B, value=x * value(B))
*(x::Real, C::Pseudoscalar) = C * x

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
*(B::AbstractScalar, C::One) = B
*(B::One, C::AbstractScalar) = C

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
*(B::AbstractScalar, C::Zero) = C
*(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
*(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) * x)
*(x::Real, B::AbstractScalar) = B * x

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
*(v::Vector{<:Real}, B::AbstractScalar) = Blade(value(B) * v)
*(B::AbstractScalar, v::Vector{<:Real}) = v * B

# ------ Specializations involving a Scalar instance

# B::Scalar, C::Scalar
*(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) * value(C))

# ------ Specializations involving a One instance

# B::One, C::One
*(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
*(B::One, C::Zero) = C
*(B::Zero, C::One) = B

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
*(B::Zero, C::Zero) = B
