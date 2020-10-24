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

export Zero

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

import Base.zero

"""
    zero(M::AbstractMultivector)
    zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the additive identity 0.
"""
zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()

# --- AbstractScalar interface functions for Zero

export value

"""
    value(B::Zero)

Return 0 (with the same precision of `B`).
"""
value(B::Zero{T}) where {T<:AbstractFloat} = T(0)

# --- Operators from the AbstractMultivector and AbstractBlade interfaces

import Base.:(+), Base.:(-)
import Base.:(*), Base.:(/)

-(B::Zero) = B
reciprocal(B::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(Inf)

+(M::AbstractMultivector, B::Zero) = M
+(B::Zero, M::AbstractMultivector) = M

-(M::AbstractMultivector, B::Zero) = M
-(B::Zero, M::AbstractMultivector) = -M

*(M::AbstractMultivector, B::Zero) = B
*(B::Zero, M::AbstractMultivector) = B

/(M::AbstractMultivector, B::Zero) = Scalar{typeof(norm(M))}(Inf)
/(B::Zero, M::AbstractMultivector) = B

# --- Operators to remove method ambiguity

# Operators between Zero instances
+(B::Zero, C::Zero) = B
-(B::Zero, C::Zero) = B
*(B::Zero, C::Zero) = B
/(B::Zero{T}, C::Zero{T}) where {T<:AbstractFloat} = Scalar{T}(NaN)
