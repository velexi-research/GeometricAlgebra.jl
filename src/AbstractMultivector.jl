#   Copyright (c) 2020-2022 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
AbstractMultivector.jl defines the AbstractMultivector type and interface
"""

# --- Exports

# Types
export AbstractMultivector

# Properties
import Base.getindex
import LinearAlgebra.norm
export dim, grades, blades, norm

# Functions
import Base.:(-), Base.reverse
export dual

import Base.:(==), Base.isapprox
import Base.iszero, Base.isone

import Base.zero, Base.one
import Base.convert

# --- Types

"""
    AbstractMultivector{<:AbstractFloat}

Supertype for all multivector types.
"""
abstract type AbstractMultivector{T<:AbstractFloat} end

# --- Public functions/methods

"""
    dim(M::AbstractMultivector)::Int

Return the dimension of the real space that `M` is embedded within.
"""
function dim end

"""
    grades(M::AbstractMultivector; collect=true)::Vector{Int}

Return the grades of the nonzero `k`-vector components of `M`. When `collect` is `false`,
an iterator over the grades is returned.
"""
function grades end

"""
    blades(M::AbstractMultivector)::Vector{<:AbstractBlade}

Return the blades that `M` is composed of.
"""
function blades end

"""
    norm(M::AbstractMultivector)::AbstractFloat

Compute the norm of `M`.
"""
function norm end

"""
    getindex(M::AbstractMultivector, k::Int)::Vector{<:AbstractBlade}

Return the `k`-vector component of `M`.
"""
function getindex end

"""
    -(M::AbstractMultivector)::AbstractMultivector

Compute the additive inverse of `M`.
"""
function -(::AbstractMultivector) end

"""
    reverse(M::AbstractMultivector)::AbstractMultivector

Compute the reverse of `M`.
"""
function reverse end

"""
    dual(M::AbstractMultivector)::AbstractMultivector

Compute the dual of `M` (relative to the unit pseudoscalar of the geometric algebra that
`M` is an element of).
"""
function dual end

# ------ Comparison methods

# B::AbstractMultivector, C::AbstractMultivector
isapprox(B::AbstractMultivector, C::AbstractMultivector) = false

# B::AbstractMultivector, x::Real
# x::Real, B::AbstractMultivector
isapprox(M::AbstractMultivector, x::Real) = false
isapprox(x::Real, M::AbstractMultivector) = false

# B::AbstractMultivector, v::Vector
# v::Vector, B::AbstractMultivector
isapprox(M::AbstractMultivector, v::Vector) = false
isapprox(v::Vector, M::AbstractMultivector) = false

# ------ Utility methods

"""
    zero(M::AbstractMultivector)::Zero
    zero(::Type{<:AbstractMultivector})::Zero
    zero(::Type{<:AbstractMultivector{T}})::Zero{T} where {T<:AbstractFloat}

Return the additive identity for the geometric algebra that `M` is an element of.
"""
zero(M::AbstractMultivector) = Zero{typeof(norm(M))}()
zero(::Type{<:AbstractMultivector}) = Zero{Float64}()
zero(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = Zero{T}()

"""
    iszero(M::AbstractMultivector)

Return `true` if `M == zero(M)`; return false otherwise.
"""
iszero(M::AbstractMultivector) = (M === zero(M))

"""
    one(M::AbstractMultivector)::One
    one(::Type{<:AbstractMultivector})::One
    one(::Type{<:AbstractMultivector{T}})::One{T} where {T<:AbstractFloat}

Return the multiplicative identity for the geometric algebra that `M` is an element of.
"""
one(M::AbstractMultivector) = One{typeof(norm(M))}()
one(::Type{<:AbstractMultivector}) = One{Float64}()
one(::Type{<:AbstractMultivector{T}}) where {T<:AbstractFloat} = One{T}()

"""
    isone(M::AbstractMultivector)

Return `true` if `M == one(M)`; return false otherwise.
"""
isone(M::AbstractMultivector) = (M === one(M))

"""
    convert(::Type{T}, M::AbstractMultivector)::AbstractMultivector{T}
        where {T<:AbstractFloat}

Convert `M` to have the floating-point precision of type `T`.
"""
function convert end

# --- Non-exported utility methods

"""
    assert_dim_equal(M::AbstractMultivector, N::AbstractMultivector)
    assert_dim_equal(M::AbstractMultivector, v::Vector)
    assert_dim_equal(v::Vector, M::AbstractMultivector)

Assert that the dimensions of `M` and `N` are the same.
"""
function assert_dim_equal(M::AbstractMultivector, N::AbstractMultivector)
    if dim(M) != dim(N)
        throw(DimensionMismatch("`dim(M)` not equal to `dim(N)`"))
    end
end

function assert_dim_equal(M::AbstractMultivector, v::Vector)
    if dim(M) != length(v)
        throw(DimensionMismatch("`dim(M)` not equal to `length(v)`"))
    end
end
assert_dim_equal(v::Vector, M::AbstractMultivector) = assert_dim_equal(M, v)
