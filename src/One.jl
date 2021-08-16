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

value(B::One{T}) where {T<:AbstractFloat} = T(1)

# --- Method definitions for AbstractBlade interface functions

reciprocal(B::One) = B
