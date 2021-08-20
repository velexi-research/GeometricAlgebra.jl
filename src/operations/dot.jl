"""
dot.jl defines methods for the dot(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

import LinearAlgebra.dot
export dot

# --- Method definitions

"""
    dot(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    ⋅(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

Compute the dot product of `M` and `N`.
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
