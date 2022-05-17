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
subtract.jl defines methods for the -(x, y) function
"""

# --- Exports

import Base.:(-)

# --- Method definitions

# M::AbstractMultivector, N::AbstractMultivector
"""
    -(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the difference between `M` and `N`.
"""
-(M::AbstractMultivector, N::AbstractMultivector) = M + -N

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
-(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) - value(C))

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
-(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) - 1)
-(B::One, C::AbstractScalar) = Scalar{typeof(value(C))}(1 - value(C))

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
-(B::AbstractScalar, C::Zero) = B
-(B::Zero, C::AbstractScalar) = -C

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
-(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) - x)
-(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x - value(B))

# ------ Specializations involving a Scalar instance

# B::Scalar, C::Scalar
-(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) - value(C))

# ------ Specializations involving a One instance

# B::One, C::One
-(B::One, C::One) = Zero{typeof(value(B))}()

# B::One, C::Zero
# B::Zero, C::One
-(B::One, C::Zero) = B
-(B::Zero, C::One) = -C

# ------ Specializations involving a Zero instance

# B::Zero, C::Zero
-(B::Zero, C::Zero) = B

# B::Zero, v::Vector
# v::Vector, B::Zero
-(B::Zero, v::Vector{<:Real}) = -v
-(v::Vector{<:Real}, B::Zero) = v
