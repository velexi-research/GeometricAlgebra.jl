"""
One.jl defines the One type and basic functions

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
export One

# Methods
import Base.one

"""
    struct One{T<:AbstractFloat} <: AbstractScalar{T}

Multiplicative identity.
"""
struct One{T<:AbstractFloat} <: AbstractScalar{T} end

"""
    One()

Alias for a One{Float64}().
"""
One() = One{Float64}()

# --- Basic functions

"""
    one(M::AbstractMultivector)
    one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()

# --- AbstractScalar interface functions for One

"""
    value(B::One)

Return 1 (with the same precision of `B`).
"""
value(B::One{T}) where {T<:AbstractFloat} = T(1)
