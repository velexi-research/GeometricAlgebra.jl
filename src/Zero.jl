"""
Zero.jl defines the Zero type and basic functions

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exported types and interface methods

# Types
export Zero

# Methods
import Base.zero

"""
    struct Zero{T<:AbstractFloat} <: AbstractScalar{T}

Additive identity.
"""
struct Zero{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    Zero()

Alias for a Zero{Float64}().
"""
Zero() = Zero{Float64}()

# --- Basic functions

"""
    zero(M::AbstractMultivector)
    zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()

# --- AbstractScalar interface functions for Zero

"""
    value(B::Zero)

Return 0 (with the same precision of `B`).
"""
value(B::Zero{T}) where {T<:AbstractFloat} = T(0)
