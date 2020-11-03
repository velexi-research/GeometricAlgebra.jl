"""
One.jl defines the One type and core methods

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
export One

# Functions
import Base.one, Base.isone

# --- Type definitions

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

# --- Method definitions for AbstractScalar interface functions

"""
    value(B::One)

Return 1 (with the same precision of `B`).
"""
value(B::One{T}) where {T<:AbstractFloat} = T(1)

# --- Method definitions for AbstractBlade interface functions

reciprocal(B::One) = B

# --- Comparison methods

isone(M::AbstractMultivector) = (M === one(M))

# --- Utility methods

"""
    one(M::AbstractMultivector)
    one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()
