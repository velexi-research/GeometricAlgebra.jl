"""
Scalar.jl defines the Scalar type and core methods

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
export Scalar

# --- Type definitions

"""
    struct Scalar{T<:AbstractFloat} <: AbstractScalar{T}

Scalar (0-blade) represented with the floating-point precision of type `T`. The
`basis` and `volume` of a Scalar are `1` and the value of the Scalar,
respectively.
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
    Scalar(value::Real)

Construct a Scalar having the specified `value`.

When the precision is not specified, the following rules are applied to set
the precision of the Scalar.

* If `value` is a floating-point value, the precision of the Scalar is inferred
  from the precision of `value`.

* If `value` is an integer, the precision of the Scalar defaults to `Float64`.
"""
Scalar(value::AbstractFloat) = Scalar{typeof(value)}(value)
Scalar(value::Integer) = Scalar(Float64(value))

# --- Method definitions for AbstractScalar interface functions

value(B::Scalar) = B.value
