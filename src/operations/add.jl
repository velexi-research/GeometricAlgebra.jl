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
add.jl defines methods for the +(x, y) function
"""

# --- Exports

import Base.:(+)

# --- Method definitions

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, N::AbstractMultivector
"""
    +(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the sum of `M` and `N`.
"""
function +(M::AbstractMultivector, N::AbstractMultivector)
    # TODO: implement this method
    return nothing
end

# M::AbstractMultivector, B::One
# B::One,M::AbstractMultivector
@inline +(M::AbstractMultivector, B::One) = Multivector(vcat([B], blades(M)))
@inline +(B::One, M::AbstractMultivector) = B + M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
@inline +(M::AbstractMultivector, B::Zero) = M
@inline +(B::Zero, M::AbstractMultivector) = M

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function +(B::Blade, C::Blade)
    assert_dim_equal(B, C)

    if grade(B) == grade(C) == 1
        return Blade(volume(B) * basis(B) + volume(C) * basis(C))
    end

    # TODO: general blade sum
end

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
@inline +(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B; value=value(B) + value(C))

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::AbstractScalar
@inline function +(B::AbstractScalar, C::AbstractScalar)
    return Scalar{typeof(value(B))}(value(B) + value(C))
end

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
@inline +(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
@inline +(B::One, C::AbstractScalar) = C + B

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
@inline +(B::AbstractScalar, C::Zero) = B
@inline +(B::Zero, C::AbstractScalar) = C

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
@inline +(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) + x)
@inline +(x::Real, B::AbstractScalar) = B + x

# ------ Specializations involving a One instance

# B::One, C::One
@inline +(B::One, C::One) = Scalar{typeof(value(B))}(2)

# B::One, C::Zero
# B::Zero, C::One
@inline +(B::One, C::Zero) = B
@inline +(B::Zero, C::One) = C

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
@inline +(B::Zero, C::Zero) = B

# B::Zero, v::Vector
# v::Vector, B::Zero
@inline +(B::Zero, v::Vector{<:Real}) = v
@inline +(v::Vector{<:Real}, B::Zero) = v
