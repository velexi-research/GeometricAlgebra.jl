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
Zero.jl defines the Zero type and core methods
"""

# --- Exports

# Types
export Zero

# --- Type definitions

"""
    struct Zero{T<:AbstractFloat} <: AbstractScalar{T}

Additive identity for a geometric algebra (extended from a real vector space of arbitrary
dimension).
"""
struct Zero{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    Zero()

Alias for a Zero{Float64}().
"""
Zero() = Zero{Float64}()

# --- Method definitions for AbstractScalar interface functions

value(B::Zero{T}) where {T<:AbstractFloat} = T(0)

# --- Method definitions for AbstractMultivector interface functions

-(B::Zero) = B

dual(B::Zero, dim::Integer) = dual_of_zero()
dual(B::Zero) = dual_of_zero()

# --- Comparison methods

# ------ ==(B, C)

import LinearAlgebra: norm

# B::Zero, v::Vector
# v::Vector, B::Zero
==(B::Zero, v::Vector{<:Real}) = (norm(v) == 0)
==(v::Vector{<:Real}, B::Zero) = (norm(v) == 0)

# B::Zero{T1}, v::Vector{T2}
isapprox(B::Zero{T1}, v::Vector{T2};
    atol::Real=0,
    rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                            T2<:AbstractFloat} =
        isapprox(value(B), norm(v), atol=atol, rtol=rtol)

# v::Vector{T1}, B::Zero{T2}
@inline isapprox(v::Vector{T1}, B::Zero{T2};
    atol::Real=0,
    rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                            T2<:AbstractFloat} =
        isapprox(B, v, atol=atol, rtol=rtol)
