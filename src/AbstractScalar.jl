"""
AbstractScalar.jl defines the AbstractScalar type and basic functions

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

export AbstractScalar

"""
    AbstractScalar{<:AbstractFloat}

Supertype for all scalar types.

Methods
-------
    value(B)::AbstractFloat
"""
abstract type AbstractScalar{T<:AbstractFloat} <: AbstractBlade{T} end

# --- AbstractMultivector interface functions for AbstractScalar

import LinearAlgebra.norm
export dim, norm, grade, basis, volume

"""
    dim(B::AbstractScalar)

Return 0.
"""
dim(B::AbstractScalar) = 0

"""
    norm(B::AbstractScalar)

Return asbolute value of `B`.
"""
norm(B::AbstractScalar) = abs(value(B))

# --- AbstractBlade interface functions for AbstractScalar

"""
    grade(B::AbstractScalar)

Return 0.
"""
grade(B::AbstractScalar) = 0

"""
    basis(B::AbstractScalar; normalized::Bool=true)

Return 1.
"""
basis(B::AbstractScalar; normalized::Bool=true) =
    normalized ? 1 : value(B)

"""
    volume(B::AbstractScalar)

Return the value of `B`.
"""
volume(B::AbstractScalar) = value(B)

"""
    sign(B::AbstractScalar)

Return the sign of the value of `B`.
"""
Base.sign(B::AbstractScalar)::Int8 = Base.sign(value(B))
