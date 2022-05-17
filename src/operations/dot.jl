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
dot.jl defines methods for the dot(x, y) function
"""

# --- Exports

import LinearAlgebra.dot
export dot

# --- Method definitions

"""
    dot(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the inner product of `M` and `N`.
"""
dot(M::AbstractMultivector, N::AbstractMultivector; left=true) =
    left ? contract_left(M, N) : nothing  # TODO

# ------ Specializations involving an AbstractMultivector instance

dot(M::AbstractMultivector, x::Real; left=true) =
    left ? contract_left(M, x) : nothing  # TODO
dot(x::Real, M::AbstractMultivector; left=true) =
    left ? contract_left(x, M) : nothing  # TODO

dot(M::AbstractMultivector, v::Vector{<:Real}; left=true) =
    left ? contract_left(M, v) : nothing  # TODO
dot(v::Vector{<:Real}, M::AbstractMultivector; left=true) =
    left ? contract_left(v, M) : nothing  # TODO
