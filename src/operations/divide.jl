#   Copyright (c) 2020-2022 Velexi Corporation
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
divide.jl defines methods for the /(x, y) function
"""

# --- Exports

import Base.:(/)

# --- Method definitions

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, B::One
# B::One, M::AbstractMultivector
/(M::AbstractMultivector, B::One) = M
/(B::One, M::AbstractMultivector) = nothing  # TODO

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
/(B::Zero, M::AbstractMultivector) = B
/(M::AbstractMultivector, B::Zero) = Scalar{typeof(norm(M))}(Inf)

# ------ Specializations involving an AbstractBlade instance

/(B::AbstractBlade, C::AbstractBlade) = B * inv(C)

# B::AbstractBlade, C::One
# B::One, C::AbstractBlade
/(B::AbstractBlade, C::One) = B
/(B::One, C::AbstractBlade) = inv(B)

# B::Zero, C::AbstractBlade
/(B::Zero, C::AbstractBlade) = B

# B::AbstractBlade, x::Real
# x::Real, B::AbstractBlade
/(B::AbstractBlade, x::Real) = Blade(B, volume=volume(B) / x)
/(x::Real, B::AbstractBlade) = x * inv(B)

# ------ Specializations involving a Blade instance

# B::Blade, C::AbstractScalar
# B::AbstractScalar, C::Blade
/(B::Blade, C::AbstractScalar) = Blade(B, volume=volume(B) / value(C))
/(B::AbstractScalar, C::Blade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=value(B) / volume(C)) :
        Blade(C, volume=-value(B) / volume(C))

# B::Blade, C::One
# B::One, C::Blade
/(B::Blade, C::One) = B
/(B::One, C::Blade) = inv(C::Blade)

# B::Blade, C::Zero
# B::Zero, C::Blade
/(B::Blade, C::Zero) = Blade(B, volume=sign(B) * Inf)
/(B::Zero, C::Blade) = B

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function /(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar{typeof(value(B))}(value(B) / value(C))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, C::Pseudoscalar
/(B::Pseudoscalar, C::AbstractScalar) = B / value(C)
/(B::AbstractScalar, C::Pseudoscalar) = value(B) / C

# B::Pseudoscalar, C::One
# B::One, C::Pseudoscalar
/(B::Pseudoscalar, C::One) = B
/(B::One, C::Pseudoscalar) = inv(C)

# B::Pseudoscalar, C::Zero
# B::Zero, C::Pseudoscalar
/(B::Pseudoscalar, C::Zero) = Pseudoscalar(B, value=sign(B) * Inf)
/(B::Zero, C::Pseudoscalar) = B

# B::Pseudoscalar, x::Real
# x::Real, B::Pseudoscalar
/(B::Pseudoscalar, x::Real) = Pseudoscalar(B, value=value(B) / x)
/(x::Real, B::Pseudoscalar) =
    mod(dim(B), 4) < 2 ?
        Pseudoscalar(B, value=x / value(B)) :
        Pseudoscalar(B, value=-x / value(B))

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
/(B::AbstractScalar, C::One) = B
/(B::One, C::AbstractScalar) = inv(C)

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
/(B::AbstractScalar, C::Zero) =
    sign(B) > 0 ?
        Scalar{typeof(norm(B))}(Inf) :
        Scalar{typeof(norm(B))}(-Inf)

/(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
/(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) / x)
/(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x / value(B))

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
/(v::Vector{<:Real}, B::AbstractScalar) = Blade(v) / value(B)

# ------ Specializations involving a Scalar instance

# B::Scalar, C::Scalar
/(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) / value(C))

# ------ Specializations involving a One instance

# B::One, C::One
/(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
/(B::One, C::Zero) = inv(C)
/(B::Zero, C::One) = B

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
/(B::Zero, C::Zero) = Scalar{typeof(value(B))}(NaN)
