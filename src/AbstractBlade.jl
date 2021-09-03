"""
AbstractBlade.jl defines the AbstractBlade type, interface, and core methods

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
export AbstractBlade

# Properties
import Base.sign
export basis, grade, volume

# Functions
export reciprocal

# --- Type definitions

"""
    AbstractBlade{<:AbstractFloat}

Supertype for all blade types.

Interface
=========

Properties
----------

    grade(B::AbstractBlade)::Int

    basis(B::AbstractBlade)::Matrix

    volume(B::AbstractBlade)::AbstractFloat

    sign(B::AbstractBlade)::Int8

Operations
----------

    dual(B::AbstractBlade, C::AbstractBlade)::AbstractBlade

    reciprocal(B::AbstractBlade)::AbstractBlade

Implementation
==============

* For the AbstractBlade type, the `volume` of a blade should encode both norm
  and orientation information.

* The return value of all methods should preserve the precision of its
  AbstractBlade arguments (when possible).

"""
abstract type AbstractBlade{T<:AbstractFloat} <: AbstractMultivector{T} end

# --- Method definitions
#
# Note: the following method definitions are no-op place holders and intended
#       to be extended.

"""
    grade(B::AbstractBlade)::Int

Return grade of `B`.
"""
function grade end

"""
    basis(B::AbstractBlade)::Matrix

Return an orthonormal basis for the subspace represented by `B`.

When `B` is an AbstractScalar, 1 (at the precision of `B`) is returned. When
`B` is a Pseudoscalar, LinearAlgebra.I is returned.
"""
function basis end

"""
    volume(B::AbstractBlade)::AbstractFloat

Return the signed volume of `B`.

When `B` is a Blade, `volume(B)` is the signed norm of the blade _relative_ to
its unit basis.
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

Notes
=====

* `dual(B, C)` is only defined if (1) `B` and `C` are extended from real
  vector spaces of the same dimension and (2) the subspace represented by `B`
  is contained in subspace represented by `C`.

* The volume of `C` is ignored.
"""
dual(B::AbstractBlade, C::AbstractBlade)::AbstractBlade =
    nothing  # no-op method to attach docstring to because an empty generic function for
             # `dual` already exists in AbstractMultivector.jl

"""
    reciprocal(B::AbstractBlade)::AbstractBlade

Compute the multiplicative inverse of blade `B`.
"""
function reciprocal end

# --- Method definitions for AbstractMultivector interface functions

grades(B::AbstractBlade) = [grade(B)]

blades(B::AbstractBlade) = Vector{AbstractBlade}([B])

norm(B::AbstractBlade) = abs(volume(B))

Base.getindex(B::AbstractBlade, k::Int) =
    k == grade(B) ? Vector{AbstractBlade}([B]) : Vector{AbstractBlade}()
