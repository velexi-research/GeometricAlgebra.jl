#   Copyright 2020 Velexi Corporation
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

`Pseudoscalar` (an ``n``-blade) represented with the floating-point precision of type `T`.
The `basis` for a `Pseudoscalar` is the standard basis for an ``n``-dimensional real vector
space. The norm and orientation of a `Pseudoscalar` are encoded in its `value`. The norm of
a `Pseudoscalar` is equal to `abs(value)` and the orientation of a `Pseudoscalar` relative
to the standard basis is equal to `sign(value)`.

Fields
------
* `dim`: the dimension of the space that the blade is embedded in

* `value`: the value of the pseudoscalar
"""
struct Pseudoscalar{T<:AbstractFloat} <: AbstractBlade{T}
    # Fields
    dim::Int
    value::T

    """
    Construct a pseudoscalar for a geometric algebra in `dim` dimensions having the
    specified `value`.
    """
    function Pseudoscalar{T}(dim::Integer, value::Real) where {T<:AbstractFloat}
        return if dim <= 0
            throw(ArgumentError("`dim` must be positive"))
        elseif value == 0
            zero(Pseudoscalar{T})
        else
            new(dim, T(value))
        end
    end
end

"""
    Pseudoscalar{T}(dim::Integer, value::Real) where {T<:AbstractFloat}

    Pseudoscalar(dim::Integer, value::Real)

Construct a pseudoscalar for a geometric algebra in `dim` dimensions having the specified
`value`.

When the precision is not specified, the following rules are applied to set the precision
of the `Pseudoscalar`.

!!! note

    * If `value` is a floating-point value, the precision of the `Pseudoscalar` is inferred
      from the precision of `value`.

    * If `value` is an integer, the precision of the Pseudoscalar is set to `Float64`.
"""
Pseudoscalar(dim::Integer, value::AbstractFloat) = Pseudoscalar{typeof(value)}(dim, value)

Pseudoscalar(dim::Integer, value::Integer) = Pseudoscalar(dim, Float64(value))

"""
    Pseudoscalar(B::Pseudoscalar{T};
                 value::Real=value(B)) where {T<:AbstractFloat}

Copy constructor. Construct a `Pseudoscalar` representing the same space as `B` having the
specified `value`.
"""
function Pseudoscalar(B::Pseudoscalar; value::Real=value(B))
    return Pseudoscalar{typeof(B.value)}(dim(B), value)
end

# --- Method definitions

"""
    value(B::Pseudoscalar)::AbstractFloat

Return the value of `B` (with the same precision as `B`).
"""
value(B::Pseudoscalar) = B.value

# --- Method definitions for AbstractBlade interface functions

import LinearAlgebra.I

function inv(B::Pseudoscalar)
    return if mod(grade(B), 4) < 2
        Pseudoscalar(B; value=1 / value(B))
    else
        Pseudoscalar(B; value=-1 / value(B))
    end
end

grade(B::Pseudoscalar) = B.dim

basis(B::Pseudoscalar) = LinearAlgebra.I

volume(B::Pseudoscalar) = value(B)

# --- Method definitions for AbstractMultivector interface functions

dim(B::Pseudoscalar) = B.dim

-(B::Pseudoscalar) = Pseudoscalar(B; value=-value(B))

Base.reverse(B::Pseudoscalar) = mod(grade(B), 4) < 2 ? B : Pseudoscalar(B; value=-value(B))

dual(B::Pseudoscalar) = Scalar(value(B))

# --- Comparison methods

==(B::Pseudoscalar, C::Pseudoscalar) = (dim(B) == dim(C)) && (value(B) == value(C))

==(B::Pseudoscalar, C::Vector) = (dim(B) == 1) && (length(C) == 1) && (value(B) == C[1])

==(B::Vector, C::Pseudoscalar) = (C == B)

function isapprox(
    B::Pseudoscalar{T1},
    C::Pseudoscalar{T2};
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : max(√eps(T1), √eps(T2)),
) where {T1<:AbstractFloat,T2<:AbstractFloat}
    return (dim(B) == dim(C)) && isapprox(value(B), value(C); atol=atol, rtol=rtol)
end

function isapprox(
    B::Pseudoscalar{T1},
    C::Vector{T2};
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : max(√eps(T1), √eps(T2)),
) where {T1<:AbstractFloat,T2<:AbstractFloat}
    return (dim(B) == 1) &&
           (length(C) == 1) &&
           isapprox(value(B), C[1]; atol=atol, rtol=rtol)
end

function isapprox(
    B::Vector{T1},
    C::Pseudoscalar{T2};
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : max(√eps(T1), √eps(T2)),
) where {T1<:AbstractFloat,T2<:AbstractFloat}
    return isapprox(C, B; atol=atol, rtol=rtol)
end

# --- Utility methods

function convert(::Type{T}, B::Pseudoscalar) where {T<:AbstractFloat}
    return T == typeof(value(B)) ? B : Pseudoscalar{T}(dim(B), value(B))
end
