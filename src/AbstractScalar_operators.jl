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

# ------ ==(B, C)

==(B::AbstractScalar, C::AbstractScalar) = (value(B) == value(C))

# Operations involving Reals
==(B::AbstractScalar, x::Real) = (x == value(B))
==(x::Real, B::AbstractScalar) = (value(B) == x)

# ------ ≈(B, C)

≈(B::AbstractScalar{T1}, C::AbstractScalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B), value(C), atol=atol, rtol=rtol)

# Operations involving Reals
≈(B::AbstractScalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)

≈(x::Real, B::AbstractScalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)

# --- Unary operators from the AbstractMultivector and AbstractBlade interfaces

# -(B)
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

# --- Binary operators from the AbstractMultivector and AbstractBlade interfaces

# ------ +(B, C)

# Operations involving One
+(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
+(B::One, C::AbstractScalar) = C + B

# Operations involving Zero
+(B::AbstractScalar, C::Zero) = B
+(B::Zero, C::AbstractScalar) = C

# Operations involving Reals
+(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) + x)
+(x::Real, B::AbstractScalar) = B + x

# ------ -(B, C)

# Operations involving One
-(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) - 1)
-(B::One, C::AbstractScalar) = -(C - B)

# Operations involving Zero
-(B::AbstractScalar, C::Zero) = B
-(B::Zero, C::AbstractScalar) = -C

# Operations involving Reals
-(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) - x)
-(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x - value(B))

# ------ *(B, C)

# Operations involving One
*(B::AbstractScalar, C::One) = B
*(B::One, C::AbstractScalar) = C

# Operations involving Zero
*(B::AbstractScalar, C::Zero) = C
*(B::Zero, C::AbstractScalar) = B

# Operations involving Reals
*(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) * x)
*(x::Real, B::AbstractScalar) = B * x

#=
# Operations involving Pseudoscalars
*(B::Scalar, C::Pseudoscalar) = Pseudoscalar(C, value=value(B) * value(C))
*(B::Pseudoscalar, C::Scalar) = C * B
=#

#=
# Operations involving Vectors
*(B::Scalar, v::Vector{<:Real}) = B * Blade(v)
*(v::Vector{<:Real}, B::Scalar) = Blade(v) * B
=#

# ------ /(B, C)

# Operations involving One
/(B::AbstractScalar, C::One) = B
/(B::One, C::AbstractScalar) = reciprocal(C)

# Operations involving Zero
/(B::AbstractScalar, C::Zero) =
    sign(B) > 0 ?
        Scalar{typeof(norm(B))}(Inf) :
        Scalar{typeof(norm(B))}(-Inf)

/(B::Zero, C::AbstractScalar) = B

# Operations involving Reals
/(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) / x)
/(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x / value(B))

# ------ wedge(B, C)

# Operations involving One
wedge(B::AbstractScalar, C::One) = B
wedge(B::One, C::AbstractScalar) = C

# Operations involving Zero
wedge(B::AbstractScalar, C::Zero) = C
wedge(B::Zero, C::AbstractScalar) = B

# Operations involving Reals
wedge(B::AbstractScalar, x::Real) = B * x
wedge(x::Real, B::AbstractScalar) = x * B

# Operations involving Vectors
wedge(B::AbstractScalar, v::Vector{<:Real}) =  Blade(value(B) * v)
wedge(v::Vector{<:Real}, B::AbstractScalar) = wedge(B, v)

# ------ contractl(B, C)

# Operations involving One
contractl(B::AbstractScalar, C::One) = B
contractl(B::One, C::AbstractScalar) = C

# Operations involving Zero
contractl(B::AbstractScalar, C::Zero) = C
contractl(B::Zero, C::AbstractScalar) = B

# Operations involving Reals
contractl(B::AbstractScalar, x::Real) = B * x
contractl(x::Real, B::AbstractScalar) = x * B

# ------ proj(B::AbstractScalar, C)

proj(B::AbstractScalar, C::AbstractScalar) = B

# C::Zero
proj(B::AbstractScalar, C::Zero) = C

# Operations involving Reals
proj(B::AbstractScalar, x::Real) = B
proj(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x)

# Operations involving Vectors
proj(v::Vector{<:Real}, B::AbstractScalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

proj(B::AbstractScalar, v::Vector{<:Real}; return_blade::Bool=true) =
    return_blade ? B : value(B)

#= REVIEW
# Operations involving Blades
proj(B::AbstractScalar, C::Blade) = B
proj(B::Blade, C::AbstractScalar) = zero(B)
=#

#= REVIEW
# Operations involving Pseudoscalars
proj(B::AbstractScalar, C::Pseudoscalar) = B
proj(B::Pseudoscalar, C::AbstractScalar) = zero(B)
=#

# ------ dual(B, C)

# B::AbstractScalar, C::AbstractScalar
dual(B::AbstractScalar, C::AbstractScalar) = B

# Operations involving Reals
dual(B::AbstractScalar, x::Real) = B
dual(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x)

# ------ dot(B, C)

# Operations involving Vectors
#=
dot(B::AbstractScalar, v::Vector{<:Real}) = B * v
dot(v::Vector{<:Real}, B::AbstractScalar) = zero(B)
=#
