"""
AbstractMultivector.jl defines the AbstractMultivector type and basic functions

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""

export AbstractMultivector

"""
    AbstractMultivector{<:AbstractFloat}

Supertype for all multivector types.

Interface
---------

Note: the return value of all methods should preserve the precision of the
AbstractMultivector instance (when possible).

### Methods

    dim(M::AbstractMultivector)::Int
    grades(M::AbstractMultivector)::Vector{Int}
    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}
    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

    norm(M::AbstractMultivector{T})::T where {T<:AbstractFloat}

### Unary Operators

    -(M::AbstractMultivector)::AbstractMultivector
    reverse(M::AbstractMultivector)::AbstractMultivector
    dual(M::AbstractMultivector)::AbstractMultivector

### Binary Operators

    +(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    -(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    *(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    /(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    wedge(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector
    ∧(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    dot(M::AbstractMultivector, N::AbstractMultivector;
        left=true)::AbstractMultivector
    ⋅(M::AbstractMultivector, N::AbstractMultivector;
        left=true)::AbstractMultivector

    contractl(M::AbstractMultivector,
              N::AbstractMultivector)::AbstractMultivector
    <(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    contractr(M::AbstractMultivector,
              N::AbstractMultivector)::AbstractMultivector
    >(M::AbstractMultivector, N::AbstractMultivector)::AbstractMultivector

    dual(M::AbstractMultivector, C::AbstractBlade)::AbstractMultivector
"""
abstract type AbstractMultivector{T<:AbstractFloat} end
