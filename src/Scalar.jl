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
Scalar.jl defines the Scalar type and core methods
"""

# --- Exports

# Types
export Scalar

# --- Type definitions

"""
    struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}

Scalar (a 0-blade) represented with the floating-point precision of type `T`. The `basis`
and `volume` of a `Scalar` are `1` and the value of the `Scalar`, respectively.
"""
struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}
    #=
        Fields
        ------
        * `value`: the value of the scalar
    =#
    value::T

    """
    Construct a Scalar having the specified `value`.
    """
    function Scalar{T}(value::Real) where {T<:AbstractFloat}
        if value == 0
            return Zero{T}()
        elseif value == 1
            return One{T}()
        end

        new(T(value))
    end
end

"""
    Scalar{T}(value::Real) where {T<:AbstractFloat}

    Scalar(value::AbstractFloat)

    Scalar(value::Integer)

Construct a scalar having the specified `value`.

When the precision is not specified, the following rules are applied to set the precision
of the `Scalar`.

!!! note

    * If `value` is a floating-point value, the precision of the `Scalar` is inferred from
      the precision of `value`.

    * If `value` is an integer, the precision of the `Scalar` is set to to `Float64`.
"""
Scalar(value::AbstractFloat) = Scalar{typeof(value)}(value)
Scalar(value::Integer) = Scalar(Float64(value))

# --- Method definitions for AbstractScalar interface functions

value(B::Scalar) = B.value
