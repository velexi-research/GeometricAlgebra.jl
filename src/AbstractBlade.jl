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
AbstractBlade.jl defines the AbstractBlade type, interface, and core methods
"""

# --- Exports

# Types
export AbstractBlade

# Properties
import Base.sign
export basis, grade, volume

# Functions
import Base.inv

# --- Types

"""
    AbstractBlade{<:AbstractFloat}

Supertype for all blade types.
"""
abstract type AbstractBlade{T<:AbstractFloat} <: AbstractMultivector{T} end

# --- Public functions/methods

"""
    grade(B::AbstractBlade)::Int

Return the grade of `B`.
"""
function grade end

"""
    basis(B::AbstractBlade)::Matrix

Return an orthonormal basis for the subspace represented by `B`.

When `B` is an `AbstractScalar`, 1 (at the precision of `B`) is returned. When `B` is a
`Pseudoscalar`, `LinearAlgebra.I` is returned.
"""
function basis end

"""
    volume(B::AbstractBlade)::AbstractFloat

Return the signed volume of `B`.

When `B` is a `Blade`, `volume(B)` is the signed norm of the blade _relative_ to its unit
basis.

!!! note

    The `volume` of a blade encodes both norm and orientation information of a blade.
"""
function volume end

"""
    sign(B::AbstractBlade)::Int8

Return the sign of `B` relative to its unit basis.
"""
sign(B::AbstractBlade)::Int8 = sign(volume(B))

"""
    dual(B::AbstractBlade, C::AbstractBlade)::AbstractBlade

Compute the dual `B` relative to the subspace represented by `C`.

!!! note

    `dual(B, C)` is only defined if

    1. `B` and `C` are extended from real vector spaces of the same dimension and

    2. the subspace represented by `B` is contained in subspace represented by `C`.

    If either of these conditions is not satisfied, an error is thrown.

!!! note

    The volume of `C` is ignored when computing `dual(B, C)`.
"""
function dual(B::AbstractBlade, C::AbstractBlade)::AbstractBlade
    # no-op method to attach docstring to because an empty generic function for `dual`
    # already exists in AbstractMultivector.jl
    return nothing
end

"""
    inv(B::AbstractBlade)::AbstractBlade

Compute the multiplicative inverse of blade `B`.
"""
function inv end

# --- Method definitions for AbstractMultivector interface functions

grades(B::AbstractBlade) = [grade(B)]

blades(B::AbstractBlade) = Vector{AbstractBlade}([B])

norm(B::AbstractBlade) = abs(volume(B))

Base.getindex(B::AbstractBlade, k::Int) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()
