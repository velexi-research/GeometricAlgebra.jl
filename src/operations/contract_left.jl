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
contract_left.jl defines methods for the contract_left(x, y) function
"""

# --- Exports

export contract_left
import Base.:(<)

# --- Operator aliases

"""
    <(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Alias for the `contract_left()` function. Compute the left contraction of `M` with `N`.
"""
<(M::AbstractMultivector, N::AbstractMultivector) = contract_left(M, N)

<(M::AbstractMultivector, x::Real) = contract_left(M, x)
<(x::Real, M::AbstractMultivector) = contract_left(x, M)

<(M::AbstractMultivector, v::Vector{<:Real}) = contract_left(M, v)
<(v::Vector{<:Real}, M::AbstractMultivector) = contract_left(v, M)

# --- Method definitions

"""
    contract_left(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the left contraction of `M` with `N`.
"""
function contract_left end

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, B::AbstractBlade
# B::AbstractBlade, M::AbstractMultivector
function contract_left(M::AbstractMultivector, B::AbstractBlade)
    return Multivector(map(C -> contract_left(C, B), blades(M)))
end

function contract_left(B::AbstractBlade, M::AbstractMultivector)
    return Multivector(map(C -> contract_left(B, C), blades(M)))
end

# M::AbstractMultivector, B::AbstractScalar
# B::AbstractScalar, M::AbstractMultivector
function contract_left(M::AbstractMultivector, B::AbstractScalar)
    return length(M[0]) > 0 ? Scalar(value(M[0][1]) * B) : zero(Scalar{typeof(norm(M))})
end

contract_left(B::AbstractScalar, M::AbstractMultivector) = B * M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
contract_left(B::Zero, M::AbstractMultivector) = zero(M)
contract_left(M::AbstractMultivector, B::Zero) = zero(M)

# M::AbstractMultivector, x::Real
# x::Real, M::AbstractMultivector
function contract_left(M::AbstractMultivector, x::Real)
    return length(M[0]) > 0 ? Scalar(value(M[0][1]) * x) : zero(Scalar{typeof(norm(M))})
end

contract_left(x::Real, M::AbstractMultivector) = x * M

# ------ Specializations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
contract_left(B::AbstractBlade, C::AbstractScalar) = zero(B)

contract_left(B::AbstractScalar, C::AbstractBlade) = Blade(C; volume=value(B) * volume(C))

# B::AbstractBlade, C::Zero
# B::Zero, C::AbstractBlade
contract_left(B::AbstractBlade, C::Zero) = zero(B)
contract_left(B::Zero, C::Blade) = zero(B)

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function contract_left(B::Blade, C::Blade)
    # --- Check arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # (B ⋅ C) = 0 if grade(B) > grade(C)
    if grade(B) > grade(C)
        return zero(B)
    end

    # --- Compute (B ⋅ C) = project(B, C) * C
    #     = volume(C) dual(project(B, C), C) * I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) project(B, C) / I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C)
    #       dual(project(B, C), C)
    #
    #     where I_C is the unit blade for the subspace represented by blade `C`
    #     that has the same orientation as basis(C).

    # Compute project(B, C) = (B ⋅ C) / C
    projection = project(B, C)

    # Compute
    #   (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) dual(project(B, C), C)
    return if mod(grade(C), 4) < 2
        volume(C) * dual(projection, C)
    else
        -volume(C) * dual(projection, C)
    end
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function contract_left(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return mod(grade(C), 4) < 2 ? dual(B, C) * volume(C) : -dual(B, C) * volume(C)
end

function contract_left(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    return zero(B)
end

# B::Blade, v::Vector
# v::Vector, B::Blade
function contract_left(B::Blade, v::Vector{<:Real})
    assert_dim_equal(v, B)
    return grade(B) > 1 ? zero(B) : Scalar(volume(B) * basis(B) ⋅ v)
end

function contract_left(v::Vector{<:Real}, B::Blade)
    # --- Check arguments

    assert_dim_equal(v, B)

    # --- Compute (v ⋅ B) = project(v, B) * B
    #     = volume(B) dual(project(v, B), B) * I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) project(v, B) / I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B)
    #       dual(project(v, B), B)
    #
    #     where I_B is the unit blade for the subspace represented by blade `B`
    #     that has the same orientation as basis(B).

    # Compute project(v, B) = (v ⋅ B) / B
    projection = project(v, B)

    # Compute
    #   (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) dual(project(v, B), B)
    return if mod(grade(B), 4) < 2
        volume(B) * dual(projection, B)
    else
        -volume(B) * dual(projection, B)
    end
end

# B::Blade, x::Real
# x::Real, B::Blade
contract_left(B::Blade, x::Real) = zero(B)
contract_left(x::Real, B::Blade) = x * B

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function contract_left(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return mod(grade(B), 4) < 2 ? Scalar(value(B) * value(C)) : Scalar(-value(B) * value(C))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, B::Pseudoscalar
contract_left(B::Pseudoscalar, C::AbstractScalar) = zero(B)
function contract_left(B::AbstractScalar, C::Pseudoscalar)
    return Pseudoscalar(C; value=value(B) * value(C))
end

# B::Pseudoscalar, C::Zero
# B::Zero, C::Pseudoscalar
contract_left(B::Pseudoscalar, C::Zero) = zero(B)
contract_left(B::Zero, C::Pseudoscalar) = zero(B)

# B::Pseudoscalar, v::Vector
# v::Vector, B::Pseudoscalar
function contract_left(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    return zero(B)
end

function contract_left(v::Vector{<:Real}, B::Pseudoscalar)
    assert_dim_equal(B, v)

    return mod(grade(B), 4) < 2 ? volume(B) * dual(Blade(v)) : -volume(B) * dual(Blade(v))
end

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
contract_left(B::AbstractScalar, C::One) = B
contract_left(B::One, C::AbstractScalar) = C

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
contract_left(B::AbstractScalar, C::Zero) = C
contract_left(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
contract_left(B::AbstractScalar, x::Real) = B * x
contract_left(x::Real, B::AbstractScalar) = x * B

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
contract_left(B::AbstractScalar, v::Vector{<:Real}) = Blade(value(B) * v)
contract_left(v::Vector{<:Real}, B::AbstractScalar) = zero(B)

# ------ Specializations involving a Scalar instance

# B::Scalar, C::Scalar
contract_left(B::Scalar, C::Scalar) = B * C

# ------ Specializations involving a One instance

# B::One, C::One
contract_left(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
contract_left(B::One, C::Zero) = C
contract_left(B::Zero, C::One) = B

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
contract_left(B::Zero, C::Zero) = B
