"""
AbstractScalar_operators.jl defines operators for the AbstractScalar type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Comparison operators

==(B::AbstractScalar, C::AbstractScalar) = (value(B) == value(C))
==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (value(B) == x)

≈(B::AbstractScalar{T1}, C::AbstractScalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B), value(C), atol=atol, rtol=rtol)

≈(B::AbstractScalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)

≈(x::Real, B::AbstractScalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)

# --- Operators from the AbstractMultivector and AbstractBlade interfaces

# ------ Unary operators

-(B::AbstractScalar) = Scalar(-value(B))

"""
    reverse(B::AbstractScalar)

Return `B` (the reverse of a scalar is itself).
"""
Base.reverse(B::AbstractScalar) = B

"""
    dual(B::AbstractScalar; dim::Integer)

Compute the dual of `B`. Note that the dimension of the embedding space must
be explicitly specified.
"""
function dual(B::AbstractScalar; dim::Union{Integer, Nothing}=nothing)
    if dim === nothing
        error("The dual of a scalar is not well-defined if `dim` is not " *
              "specified")
    end

    mod(dim, 4) < 2 ?
        Pseudoscalar(dim, value(B)) :
        Pseudoscalar(dim, -value(B))
end

reciprocal(B::AbstractScalar) = 1 / B

# --- Special cases

# Operations involving AbstractMultivectors
*(B::Scalar, M::AbstractMultivector) = value(B) * M
*(M::AbstractMultivector, B::Scalar) = B * M

contractl(B::Scalar, M::AbstractMultivector) = B * M
contractl(M::AbstractMultivector, B::Scalar) =
    length(M[0]) > 0 ?
        Scalar(value(M[0][1]) * B) :
        zero(Scalar{typeof(norm(M))})

# Operations involving Reals
+(B::AbstractScalar, C::Real) = Scalar{typeof(value(B))}(value(B) + C)
+(B::Real, C::AbstractScalar) = C + B

-(B::AbstractScalar, C::Real) = Scalar{typeof(value(B))}(value(B) - C)
-(B::Real, C::AbstractScalar) = Scalar{typeof(value(C))}(B - value(C))

*(B::AbstractScalar, C::Real) = Scalar{typeof(value(B))}(value(B) * C)
*(B::Real, C::AbstractScalar) = C * B

/(B::AbstractScalar, C::Real) = Scalar{typeof(value(B))}(value(B) / C)
/(B::Real, C::AbstractScalar) = Scalar{typeof(value(C))}(B / value(C))

wedge(B::AbstractScalar, C::Real) = B * C
wedge(B::Real, C::AbstractScalar) = B * C

contractl(B::AbstractScalar, C::Real) = B * C
contractl(B::Real, C::AbstractScalar) = B * C

proj(B::AbstractScalar, C::Real) = B
proj(B::Real, C::AbstractScalar) = B

dual(B::AbstractScalar, C::Real) = B
dual(B::Real, C::AbstractScalar) = B
