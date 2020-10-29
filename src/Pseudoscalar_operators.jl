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

*(B::Pseudoscalar, C::Pseudoscalar) = dot(B, C)

#= REVIEW
# Operations involving Blades
*(B::Blade, C::Pseudoscalar) = B ⋅ C
*(B::Pseudoscalar, C::Blade) = zero(B)

# Operations involving Vectors
*(B::Pseudoscalar, v::Vector{<:Real}) = zero(B)
*(v::Vector{<:Real}, B::Pseudoscalar) = Blade(v) * B

# Operations involving Reals
*(x::Real, B::Pseudoscalar) = Pseudoscalar(B, value=x * value(B))
*(B::Pseudoscalar, x::Real) = x * B
=#

# ------ /(B, C)

/(B::Pseudoscalar, C::Pseudoscalar) =
    Scalar{typeof(value(B))}(value(B) / value(C))

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

#= REVIEW
# Operations involving Blades
function contractl(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function contractl(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)

    mod(grade(C), 4) < 2 ?
        dual(B, C) * volume(C) :
       -dual(B, C) * volume(C)
end

# Operations involving Vectors
function contractl(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(v, B)
    zero(B)
end

function contractl(v::Vector{<:Real}, B::Pseudoscalar)
    assert_dim_equal(v, B)
    mod(grade(B), 4) < 2 ?
        value(B) * dual(Blade(v)) :
       -value(B) * dual(Blade(v))
end

# Operations involving Reals
contractl(x::Real, B::Pseudoscalar) = x * B
contractl(B::Pseudoscalar, x::Real) = zero(B)
=#

# ------ proj(B, C)

function proj(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

#= REVIEW
# Operations involving Blades
function proj(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function proj(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end
=#

# ------ dual(B, C)

function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar{typeof(value(B))}(value(B))
end

#= REVIEW
# Operations involving Blades
function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function dual(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    dual(B)
end
=#
