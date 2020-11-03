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

# --- Method definitions

"""
    dot(M, N)
    M ⋅ N

Compute the dot product of the multivector `M` with the multivector `N`.
"""
dot(M::AbstractMultivector, N::AbstractMultivector; left=True) =
    left ? contractl(M, N) : nothing  # TODO

dot(M::AbstractMultivector, x::Real; left=True) =
    left ? contractl(M, x) : nothing  # TODO
dot(x::Real, M::AbstractMultivector; left=True) =
    left ? contractl(x, M) : nothing  # TODO

dot(M::AbstractMultivector, v::Vector{<:Real}; left=True) =
    left ? contractl(M, v) : nothing  # TODO
dot(v::Vector{<:Real}, M::AbstractMultivector; left=True) =
    left ? contractl(v, M) : nothing  # TODO
