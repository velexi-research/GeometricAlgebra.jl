"""
Zero.jl defines the Zero type and core methods

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

# Types
export Zero

# Functions
import Base.zero, Base.iszero

# --- Type definitions

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

# --- Method definitions for AbstractScalar interface functions

"""
    value(B::Zero)

Return 0 (with the same precision of `B`).
"""
value(B::Zero{T}) where {T<:AbstractFloat} = T(0)

# --- Method definitions for AbstractMultivector interface functions

-(B::Zero) = B

dual(B::Zero; dim::Union{Integer, Nothing}=nothing) = dual_of_zero()

# --- Comparison methods

iszero(M::AbstractMultivector) = (M === zero(M))

# --- Utility methods

"""
    zero(M::AbstractMultivector)
    zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()
