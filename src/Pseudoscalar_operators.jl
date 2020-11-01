"""
Pseudoscalar_operators.jl defines the Pseudoscalar type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Comparison operators

==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

≈(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && ≈(value(B), value(C), atol=atol, rtol=rtol)

# --- Unary operators from the AbstractMultivector and AbstractBlade interfaces

-(B::Pseudoscalar) = Pseudoscalar(B, value=-value(B))

reciprocal(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Pseudoscalar(B, value=1 / value(B)) :
        Pseudoscalar(B, value=-1 / value(B))

Base.reverse(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?  B : Pseudoscalar(B, value=-value(B))

dual(B::Pseudoscalar) = Scalar(value(B))

# --- Binary operators from the AbstractMultivector and AbstractBlade interfaces

# ------ +(B, C)

+(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) + value(C))

# ------ -(B, C)

-(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) - value(C))

# ------ *(B, C)

*(B::Pseudoscalar, C::Pseudoscalar) = contractl(B, C)

# Operations involving AbstractScalars
*(B::Pseudoscalar, C::AbstractScalar) =
    Pseudoscalar(B, value=value(B) * value(C))

*(B::AbstractScalar, C::Pseudoscalar) = C * B

# Operations involving One
*(B::Pseudoscalar, C::One) = B
*(B::One, C::Pseudoscalar) = C

# Operations involving Zero
*(B::Pseudoscalar, C::Zero) = C
*(B::Zero, C::Pseudoscalar) = B

# Operations involving Reals
*(B::Pseudoscalar, x::Real) = Pseudoscalar(B, value=x * value(B))
*(x::Real, C::Pseudoscalar) = C * x

# ------ /(B, C)

function /(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar{typeof(value(B))}(value(B) / value(C))
end

# Operations involving AbstractScalar
/(B::Pseudoscalar, C::AbstractScalar) = B / value(C)
/(B::AbstractScalar, C::Pseudoscalar) = value(B) / C

# Operations involving One
/(B::Pseudoscalar, C::One) = B
/(B::One, C::Pseudoscalar) = 1 / C

# Operations involving Zero
/(B::Pseudoscalar, C::Zero) = Pseudoscalar(B, value=sign(B) * Inf)
/(B::Zero, C::Pseudoscalar) = B

# Operations involving Reals
/(B::Pseudoscalar, x::Real) = Pseudoscalar(B, value=value(B) / x)
/(x::Real, B::Pseudoscalar) =
    mod(dim(B), 4) < 2 ?
        Pseudoscalar(B, value=x / value(B)) :
        Pseudoscalar(B, value=-x / value(B))

# ------ wedge(B, C)

function wedge(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    zero(B)
end

# Operations involving Vectors
function wedge(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    zero(B)
end

wedge(v::Vector{<:Real}, B::Pseudoscalar) = B ∧ v

# ------ contractl(B, C)

function contractl(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))
end

# ------ proj(B, C)

function proj(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

# Operations involving AbstractScalars
proj(B::Pseudoscalar, C::AbstractScalar; return_blade::Bool=true) = zero(B)

# ------ dual(B, C)

# B::Pseudoscalar, C::Pseudoscalar
function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar{typeof(value(B))}(value(B))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, C::Pseudoscalar
dual(B::Pseudoscalar, C::AbstractScalar) = zero(B)

dual(B::AbstractScalar, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=value(B)) :
        Pseudoscalar(C, value=-value(B))

# Operations involving Reals
dual(B::Pseudoscalar, x::Real) = zero(B)
dual(x::Real, B::Pseudoscalar) = Scalar{typeof(value(B))}(x)
