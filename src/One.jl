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
One.jl defines the One type and core methods
"""

# --- Exports

# Types
export One

# --- Type definitions

"""
    struct One{T<:AbstractFloat} <: AbstractScalar{T}

Multiplicative identity for a geometric algebra (extended from a real vector
space of arbitrary dimension).
"""
struct One{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    One()

Alias for a One{Float64}().
"""
One() = One{Float64}()

# --- Method definitions for AbstractScalar interface functions

value(B::One{T}) where {T<:AbstractFloat} = T(1)

# --- Method definitions for AbstractBlade interface functions

reciprocal(B::One) = B
