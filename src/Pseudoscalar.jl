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
Pseudoscalar.jl defines the Pseudoscalar type and core methods
"""

# --- Exports

# Types
export Pseudoscalar

# --- Type definitions

"""
    struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade

Pseudoscalar (an ``(n-1)``-blade) represented with the floating-point precision
of type `T`. The `basis` for a Pseudoscalar is the standard basis for ``R^n``.
The norm and orientation of a Pseudoscalar are encoded in its `value`. The
norm of a Pseudoscalar is equal to `abs(value)` and the orientation of a
Pseudoscalar relative to the standard basis for ``R^n`` is equal to
`sign(value)`.
"""
struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade{T}
    #=
      Fields
      ------
      * `dim`: the dimension of the space that the blade is embedded in

      * `value`: the value of the pseudoscalar
    =#
    dim::Int
    value::T

    """
    Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having
    the specified `value`.
    """
    Pseudoscalar{T}(dim::Integer, value::Real) where {T<:AbstractFloat} =
        dim <= 0 ?
            throw(ArgumentError("`dim` must be positive")) :
            value == 0 ? zero(Pseudoscalar{T}) : new(dim, T(value))
end

"""
    Pseudoscalar(dim::Integer, value::Real)

Construct a Pseudoscalar for a geometric algebra in `dim` dimensions having the
specified `value`.

When the precision is not specified, the following rules are applied to set
the precision of the Pseudoscalar.

* If `value` is a floating-point value, the precision of the Pseudoscalar is
  inferred from the precision of `value`.

* If `value` is an integer, the precision of the Pseudoscalar defaults to
  `Float64`.
"""
Pseudoscalar(dim::Integer, value::AbstractFloat) =
    Pseudoscalar{typeof(value)}(dim, value)

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))

"""
    Pseudoscalar(B::Pseudoscalar{T};
                 value::Real=value(B)) where {T<:AbstractFloat}

Copy constructor. Construct a Pseudoscalar representing the same space as
`B` having the specified `value`.
"""
Pseudoscalar(B::Pseudoscalar; value::Real=value(B)) =
    Pseudoscalar{typeof(B.value)}(dim(B), value)

# --- Method definitions

"""
    value(B::Pseudoscalar)::AbstractFloat

Return the value of `B` (with the same precision as `B`).
"""
value(B::Pseudoscalar) = B.value

# --- Method definitions for AbstractBlade interface functions

import LinearAlgebra.I

reciprocal(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Pseudoscalar(B, value=1 / value(B)) :
        Pseudoscalar(B, value=-1 / value(B))

grade(B::Pseudoscalar) = B.dim

basis(B::Pseudoscalar) = LinearAlgebra.I

volume(B::Pseudoscalar) = value(B)

# --- Method definitions for AbstractMultivector interface functions

dim(B::Pseudoscalar) = B.dim

-(B::Pseudoscalar) = Pseudoscalar(B, value=-value(B))

Base.reverse(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?  B : Pseudoscalar(B, value=-value(B))

dual(B::Pseudoscalar) = Scalar(value(B))

# --- Comparison methods

==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

==(B::Pseudoscalar, C::Vector) = 
    (dim(B) == 1) && (length(C) == 1) && (value(B) == C[1])

==(B::Vector, C::Pseudoscalar) = (C == B)

isapprox(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && isapprox(value(B), value(C), atol=atol, rtol=rtol)

isapprox(B::Pseudoscalar{T1}, C::Vector{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == 1) && (length(C) == 1) && isapprox(value(B), C[1], atol=atol, rtol=rtol)

isapprox(B::Vector{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    isapprox(C, B, atol=atol, rtol=rtol)

# --- Utility methods

convert(::Type{T}, B::Pseudoscalar) where {T<:AbstractFloat} =
    T == typeof(value(B)) ? B : Pseudoscalar{T}(dim(B), value(B))
