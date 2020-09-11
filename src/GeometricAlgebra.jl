"""
The GeometricAlgebra.jl module defines geometric algebra types and functions.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
module GeometricAlgebra

# --- Imports

import LinearAlgebra


# --- Submodules

include("blade.jl")


# --- Functions

# Exports
export zero, one


"""
    zero(B::AbstractBlade{T}) where {T<:AbstractFloat}
    zero(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}}
    zero(::Type{Blade})
    zero(::Type{Scalar})

Return Zero (the additive identity). When the precision is not explicitly
specified, it defaults to Float64.
"""
zero(B::AbstractBlade{T}) where {T<:AbstractFloat} = Zero{T}()
zero(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}} = Zero{T}()
zero(::Type{Blade}) = Zero{Float64}()
zero(::Type{Scalar}) = Zero{Float64}()


"""
    one(B::AbstractBlade{T}) where {T<:AbstractFloat}
    one(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}}
    one(::Type{Blade})
    one(::Type{Scalar})

Return One (the multiplicative identity). When the precision is not explicitly
specified, it defaults to Float64.
"""
one(B::AbstractBlade{T}) where {T<:AbstractFloat} = One{T}()
one(::Type{S}) where {T<:AbstractFloat, S<:AbstractBlade{T}} = One{T}()
one(::Type{Blade}) = One{Float64}()
one(::Type{Scalar}) = One{Float64}()

end  # End of GeometricAlgebra.jl module
