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

export One

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

import Base.one

"""
    one(M::AbstractMultivector)
    one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat}

Return the multiplicative identity 1.
"""
one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()

# --- AbstractScalar interface functions for One

export value

"""
    value(B::One)

Return 1 (with the same precision of `B`).
"""
value(B::One{T}) where {T<:AbstractFloat} = T(1) 

# --- Operators from the AbstractMultivector and AbstractBlade interfaces

import Base.:(+), Base.:(-)
import Base.:(*), Base.:(/)

-(B::One{T}) where {T<:AbstractFloat} = Scalar{T}(-1)
reciprocal(B::One) = B

+(B::One, M::AbstractMultivector) = Multivector(vcat([B], blades(M)))
+(M::AbstractMultivector, B::One) = B + M

-(M::AbstractMultivector, B::One) = Multivector(vcat([-B], blades(M)))
-(B::One, M::AbstractMultivector) = -(M - B)

*(M::AbstractMultivector, B::One) = M
*(B::One, M::AbstractMultivector) = M

/(M::AbstractMultivector, B::One) = M
/(B::One, M::AbstractMultivector) = nothing  # TODO

# --- Operators to remove method ambiguity

# Operators between One instances
+(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Scalar{T}(2)
-(B::One{T}, C::One{T}) where {T<:AbstractFloat} = Zero{T}()
*(B::One, C::One) = B
/(B::One, C::One) = B

# Operators between One and Zero instances
+(B::Zero, C::One) = C
+(B::One, C::Zero) = B

-(B::Zero, C::One) = -C
-(B::One, C::Zero) = -B

*(B::Zero, C::One) = B
*(B::One, C::Zero) = C

/(B::Zero, C::One) = B
/(B::One, C::Zero) = reciprocal(C)

# Operators between One and AbstractScalar instances
+(B::AbstractScalar, C::One) = Scalar(B, value=value(B) + 1)
+(B::One, C::AbstractScalar) = C + B

-(B::AbstractScalar, C::One) = B + -C
-(B::One, C::AbstractScalar) = B + -C
