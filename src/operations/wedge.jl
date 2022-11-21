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
wedge.jl defines methods for the wedge(x, y) function
"""

# --- Exports

export wedge, ∧

# --- Operator aliases

"""
    ∧(M, N)::AbstractMultivector

Alias for the [`wedge()`](@ref wedge) function. Compute the exterior product of `M` and `N`.
"""
∧(M, N) = wedge(M, N)

# --- Method definitions

"""
    wedge(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    wedge(M::AbstractMultivector, N::Vector{<:Real})::AbstractMultivector
    wedge(M::Vector{<:Real}, N::AbstractMultivector)::AbstractMultivector
    wedge(M::Vector{<:Real}, N::Vector{<:Real})::AbstractMultivector

    wedge(M::AbstractMultivector, N::Real)::AbstractMultivector
    wedge(M::Real, N::AbstractMultivector)::AbstractMultivector

Compute the exterior product of `M` and `N`.
"""
function wedge end

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, B::AbstractBlade
# B::AbstractBlade, M::AbstractMultivector
function wedge(M::AbstractMultivector, B::AbstractBlade)
    return Multivector(map(C -> wedge(C, B), blades(M)))
end

function wedge(B::AbstractBlade, M::AbstractMultivector)
    return Multivector(map(C -> wedge(B, C), blades(M)))
end

# M::AbstractMultivector, v::Vector
# v::Vector, M::AbstractMultivector
function wedge(M::AbstractMultivector, v::Vector{<:Real})
    return Multivector(map(C -> wedge(C, v), blades(M)))
end

function wedge(v::Vector{<:Real}, M::AbstractMultivector)
    return Multivector(map(C -> wedge(v, C), blades(M)))
end

# M::AbstractMultivector, B::AbstractScalar
# B::AbstractScalar, M::AbstractMultivector
wedge(B::AbstractMultivector, C::AbstractScalar) = B * C
wedge(B::AbstractScalar, C::AbstractMultivector) = B * C

# M::AbstractMultivector, B::One
# B::One, M::AbstractMultivector
wedge(M::AbstractMultivector, B::One) = M
wedge(B::One, M::AbstractMultivector) = M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
wedge(B::Zero, M::AbstractMultivector) = B
wedge(M::AbstractMultivector, B::Zero) = B

# M::AbstractMultivector, x::Real
# x::Real, M::AbstractMultivector
wedge(M::AbstractMultivector, x::Real) = x * M
wedge(x::Real, M::AbstractMultivector) = x * M

# ------ Specializations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
wedge(B::AbstractBlade, C::AbstractScalar) = B * C
wedge(B::AbstractScalar, C::AbstractBlade) = B * C

# B::AbstractBlade, C::One
# B::One, C::AbstractBlade
wedge(B::AbstractBlade, C::One) = B
wedge(B::One, C::AbstractBlade) = C

# B::AbstractBlade, C::Zero
# B::Zero, C::AbstractBlade
wedge(B::AbstractBlade, C::Zero) = zero(B)
wedge(B::Zero, C::AbstractBlade) = zero(C)

# B::AbstractBlade, x::Real
# x::Real, B::AbstractBlade
wedge(B::AbstractBlade, x::Real) = B * x
wedge(x::Real, B::AbstractBlade) = x * B

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function wedge(B::Blade, C::Blade)
    # --- Check arguments

    assert_dim_equal(B, C)

    if grade(B) + grade(C) > dim(B)
        return zero(B)
    end

    # --- Construct new blade

    return volume(B) * volume(C) * Blade(hcat(basis(B), basis(C)))
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function wedge(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return zero(B)
end

wedge(B::Pseudoscalar, C::Blade) = wedge(C, B)

# B::Blade, v::Vector
# v::Vector, B::Blade
function wedge(B::Blade, v::Vector{<:Real})
    assert_dim_equal(B, v)

    # Note: volume(B) is incorporated into the norm of `v`
    return Blade(hcat(basis(B), volume(B) * v))
end

function wedge(v::Vector{<:Real}, B::Blade)
    assert_dim_equal(v, B)

    # Note: volume(B) is incorporated into the norm of `v`
    return Blade(hcat(volume(B) * v, basis(B)))
end

# v::Vector, w::Vector
wedge(v::Vector{<:Real}, w::Vector{<:Real}) = Blade(hcat(v, w))

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function wedge(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return zero(B)
end

# B::Pseudoscalar, v::Vector
# v::Vector, B::Pseudoscalar
function wedge(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    return zero(B)
end

wedge(v::Vector{<:Real}, B::Pseudoscalar) = B ∧ v

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
wedge(B::AbstractScalar, C::One) = B
wedge(B::One, C::AbstractScalar) = C

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
wedge(B::AbstractScalar, C::Zero) = C
wedge(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
wedge(B::AbstractScalar, x::Real) = B * x
wedge(x::Real, B::AbstractScalar) = x * B

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
wedge(B::AbstractScalar, v::Vector{<:Real}) = value(B) * Blade(v)
wedge(v::Vector{<:Real}, B::AbstractScalar) = wedge(B, v)

# ------ Specializations involving a Scalar instance

# B::Scalar, C::Scalar
wedge(B::Scalar, C::Scalar) = B * C

# ------ Specializations involving a One instance

# B::One, C::One
wedge(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
wedge(B::One, C::Zero) = C
wedge(B::Zero, C::One) = B

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
wedge(B::Zero, C::Zero) = B
