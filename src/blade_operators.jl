"""
The blade_operators.jl submodule defines operations on subtypes of
AbstractBlade.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
import LinearAlgebra


# --- Comparison operators

# Imports
import Base.:(==), Base.:(≈)

"""
    ==(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are equal; otherwise, return false.
"""
# Scalar comparison
==(B::Scalar, C::Scalar) = (value(B) == value(C))
==(B::Scalar, x::Real) = (x == value(B))
==(x::Real, B::Scalar) = (B == x)

# Blade comparison
==(B::Blade, C::Blade) =
    dim(B) == dim(C) && grade(B) == grade(C) &&
    volume(B) == volume(C) && basis(B) == basis(C)

# Pseudoscalar comparison
==(B::Pseudoscalar, C::Pseudoscalar) =
    (dim(B) == dim(C)) && (value(B) == value(C))

"""
    ≈(B::AbstractBlade, C::AbstractBlade)

Return true if B and C are approximatly equal; otherwise, return false.
"""
# Scalar comparison
≈(B::Scalar{T1}, C::Scalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    ≈(value(B), value(C), atol=atol, rtol=rtol)

≈(B::Scalar{T}, x::Real;
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(x, value(B), atol=atol, rtol=rtol)

≈(x::Real, B::Scalar{T};
  atol::Real=0, rtol::Real=atol>0 ? 0 : sqrt(eps(T))) where {T<:AbstractFloat} =
    ≈(B, x, atol=atol, rtol=rtol)

# Blade comparison
function ≈(B::Blade{T1}, C::Blade{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat}
    # Check dim, grade, and norm are equal
    if dim(B) != dim(C) || grade(B) != grade(C) ||
        ≉(norm(B), norm(C), atol=atol, rtol=rtol)

        return false
    end

    # Check that B and C represent the same space
    projection = LinearAlgebra.det(transpose(basis(B)) * basis(C))
    if ≉(abs(projection), 1, atol=atol, rtol=rtol)
        return false
    end

    # Check that B and C have the same orientation
    return sign(B) * sign(C) == sign(projection)
end

# Pseudoscalar comparison
≈(B::Pseudoscalar{T1}, C::Pseudoscalar{T2};
  atol::Real=0,
  rtol::Real=atol>0 ? 0 : max(√eps(T1), √eps(T2))) where {T1<:AbstractFloat,
                                                          T2<:AbstractFloat} =
    (dim(B) == dim(C)) && ≈(value(B), value(C), atol=atol, rtol=rtol)

# Default AbstractBlades comparisons
≈(B::AbstractBlade, C::AbstractBlade) = false
≈(B::AbstractBlade, C::Real) = false
≈(B::Real, C::AbstractBlade) = false


# --- Core Blade operations

import Base.:(*)
export ∧, outer
export project, dual

"""
    *(B::Union{AbstractBlade}, C::Union{Scalar, Real})
    *(B::Union{Scalar, Real}, C::Union{AbstractBlade})

Return the product of `B` and `C`.
"""
*(x::Real, B::Scalar) = Scalar(B, value=x * value(B))
*(B::Scalar, x::Real) = x * B
*(B::Scalar, C::Scalar) = Scalar(B, value=value(B) * value(C))

*(x::Real, B::Blade) = Blade(B, volume=x * volume(B))
*(B::Blade, x::Real) = x * B
*(B::Scalar, C::Blade) = value(B) * C
*(B::Blade, C::Scalar) = C * B

*(x::Real, B::Pseudoscalar) = Pseudoscalar(B, value=x * value(B))
*(B::Pseudoscalar, x::Real) = x * B
*(B::Scalar, C::Pseudoscalar) = Pseudoscalar(C, value=value(B) * value(C))
*(B::Pseudoscalar, C::Scalar) = C * B

"""
    ∧(B::Union{AbstractBlade, Vector, Real},
      C::Union{AbstractBlade, Vector, Real})

    outer(B::Union{AbstractBlade, Vector, Real},
          C::Union{AbstractBlade, Vector, Real})

Return the outer product of `B` and `C`.
"""
# Outer products involving Scalars
∧(x::Real, B::Scalar) = x * B
∧(B::Scalar, x::Real) = B * x
∧(B::Scalar, C::Scalar) = B * C

∧(x::Real, B::Blade) = x * B
∧(B::Blade, x::Real) = B * x
∧(B::Scalar, C::Blade) = B * C
∧(B::Blade, C::Scalar) = B * C

∧(x::Real, B::Pseudoscalar) = x * B
∧(B::Pseudoscalar, x::Real) = B * x
∧(B::Scalar, C::Pseudoscalar) = B * C
∧(B::Pseudoscalar, C::Scalar) = B * C

# Outer products involving Pseudoscalars
function ∧(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    zero(B)
end

function ∧(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

∧(B::Blade, C::Pseudoscalar) = C ∧ B

# Outer product between Blades
function ∧(B::Blade, C::Blade)
    assert_dim_equal(B, C)
    Blade(hcat(basis(B), basis(C)), volume=volume(B) * volume(C))
end

# Outer product between vectors and Blades
∧(v::Vector{<:Real}, B::Blade) = Blade(v) ∧ B
∧(B::Blade, v::Vector{<:Real}) = B ∧ Blade(v)

# Outer product between vectors
∧(v::Vector{<:Real}, w::Vector{<:Real}) = Blade(hcat(v, w))

# Function aliases
outer(x::AbstractBlade, y::AbstractBlade) = x ∧ y
outer(x::AbstractBlade, y::Union{Vector{<:Real}, Real}) = x ∧ y
outer(x::Union{Vector{<:Real}, Real}, y::AbstractBlade) = x ∧ y

"""
    project(v::Vector, B::AbstractBlade; return_blade::Bool=true)
    project(B::AbstractBlade, v::Vector; return_blade::Bool=true)

Return the projection of vector `v` onto the subspace represented by blade `B`.

When `return_blade` is true, the return value is an AbstractBlade. Otherwise,
the return value is a Real (if the result is a scalar) or a Vector.
"""
project(v::Vector{<:Real}, B::Scalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

project(B::Scalar, v::Vector{<:Real}; return_blade::Bool=true) =
    return_blade ? B : value(B)

project(v::Vector{<:Real}, B::Pseudoscalar; return_blade::Bool=true) =
    length(v) != dim(B) ?
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`")) :
        return_blade ? Blade(v) : v

project(B::Pseudoscalar, v::Vector{<:Real}; return_blade::Bool=true) =
    length(v) != dim(B) ?
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`")) :
        return_blade ? zero(B) : 0

function project(v::Vector{<:Real}, B::Blade; return_blade::Bool=true)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    projection = (grade(B) == 1) ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) :
        basis(B) * transpose(transpose(v) * basis(B))

    return_blade ? Blade(projection) : projection
end

function project(B::Blade, v::Vector{<:Real}; return_blade::Bool=true)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    projection = (grade(B) == 1) ?
        v * LinearAlgebra.dot(v, basis(B)) : 0

    return_blade ? Blade(projection) : projection
end

"""
    project(B::AbstractBlade, C::AbstractBlade)

Return the projection of `B` onto the subspace represented by blade `C`.
"""
# Projections involving Scalars
project(B::Scalar, C::Scalar) = B
project(B::Scalar, C::Blade) = B
project(B::Blade, C::Scalar) = zero(B)
project(B::Scalar, C::Pseudoscalar) = B
project(B::Pseudoscalar, C::Scalar) = zero(B)

# Projections involving Pseudoscalars
function project(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

function project(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function project(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

# Projection of one Blade onto another Blade
function project(B::Blade, C::Blade)
    # --- Validate arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # When grade(B) > grade(C), the projection is zero.
    if grade(B) > grade(C)
        return zero(B)
    end

    # Compute the projections of basis(B) onto basis(C)
    projections = Matrix{typeof(volume(B))}(undef, dim(B), grade(B))
    for i in 1:grade(B)
        projections[:, i] = project(basis(B)[:, i], C, return_blade=false)
    end

    # Encode volume(B) in the first projection vector
    projections[:, 1] *= volume(B)

    # Compute projection (using fact that projection is an outermorphism)
    Blade(projections)
end

"""
    dual(B::Union{Blade, Pseudoscalar})
    dual(B::Scalar, dim::Integer)

Return the dual `B` (relative to the space that the geometric algebra is
extended from). Note that when `B` is a Scalar, the dimension of the embedding
space must be explicitly specified.
"""
dual(B::Scalar, dim::Integer) =
    mod(dim, 4) < 2 ?
        Pseudoscalar(dim, value(B)) :
        Pseudoscalar(dim, -value(B))

dual(B::Pseudoscalar) = Scalar(value(B))

function dual(B::Blade)
    # --- Extend basis(B) to an orthonormal basis for entire space.

    F = LinearAlgebra.qr(basis(B))

    # --- Compute volume of dual

    # Account for orientation of Q relative to orientation of I formed from
    # standard basis
    dual_volume = volume(B) * sign(LinearAlgebra.det(F.Q))

    # Account for orientation of first grade(B) columns of Q relative to
    # orientation of basis(B)
    if prod(LinearAlgebra.diag(F.R)) < 0
        dual_volume = -dual_volume
    end

    # Account for sign of I^{-1} relative to I
    if mod(dim(B), 4) >= 2
        dual_volume = -dual_volume
    end

    # Account for reversals required to eliminate B
    if mod(grade(B), 4) >= 2
        dual_volume = -dual_volume
    end

    Blade{typeof(volume(B))}(dim(B), dim(B) - grade(B),
                             F.Q[:, grade(B) + 1:end], dual_volume,
                             enforce_constraints=false,
                             copy_basis=false)
end

"""
    dual(B::AbstractBlade, C::AbstractBlade)

Return the dual `B` relative to the subspace represented by `C`.

Notes
-----
* `dual(B, C)` is only defined if (1) `B` and `C` are extended from real
  vector spaces of the same dimension and (2) the subspace represented by `B`
  is contained in subspace represented by `C`.

* The volume of `C` is ignored.
"""
# Duals involving Scalars
dual(B::Scalar, C::Scalar) = B

dual(B::Scalar, C::Blade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=value(B)) :
        Blade(C, volume=-value(B))

dual(B::Blade, C::Scalar) = zero(B)

dual(B::Scalar, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=value(B)) :
        Pseudoscalar(C, value=-value(B))

dual(B::Pseudoscalar, C::Scalar) = zero(B)

# Duals involving Pseudoscalars
function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar(value(B))
end

function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function dual(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    dual(B)
end

# Duals involving Blades
function dual(B::Blade, C::Blade)
    # --- Validate arguments

    # Check that B and C are extended from the real vector spaces of the same
    # dimension
    assert_dim_equal(B, C)

    # Check that B is contained in C
    projection_coefficients = transpose(basis(C)) * basis(B)
    if LinearAlgebra.norm(projection_coefficients)^2 ≉ grade(B)
        error("`B` not contained in `C`")
    end

    # --- Handle edge cases

    # Subspaces represented by B and C are the same
    if grade(B) == grade(C)
        dual_sign = mod(grade(B), 4) < 2 ? 1 : -1

        dual_volume = LinearAlgebra.det(projection_coefficients) > 0 ?
            dual_sign * volume(B) : -dual_sign * volume(B)

        return Scalar{typeof(volume(B))}(dual_volume)
    end

    # --- Extend basis(B) to an orthonormal basis for entire subspace
    #     represented by basis(C)

    F = LinearAlgebra.qr(projection_coefficients)

    # --- Compute volume of dual

    # Account for orientation of Q relative to orientation of basis(C)
    dual_volume = volume(B) * sign(LinearAlgebra.det(F.Q))

    # Account for orientation of first grade(B) columns of Q relative to
    # orientation of basis(B)
    if prod(LinearAlgebra.diag(F.R)) < 0
        dual_volume = -dual_volume
    end

    # Account for sign of I_C^{-1} relative to I_C
    if mod(grade(C), 4) >= 2
        dual_volume = -dual_volume
    end

    # Account for reversals required to eliminate B
    if mod(grade(B), 4) >= 2
        dual_volume = -dual_volume
    end

    # --- Construct Blade in embedding space

    Blade(basis(C) * F.Q[:, grade(B) + 1:end], volume=dual_volume)
end


# --- Unary operations

# Imports
import Base.:(-)

# Exports
export reciprocal, reverse

"""
    -(B::AbstractBlade)

Return the additive inverse of `B`.
"""
-(B::Scalar) = Scalar(B, value=-value(B))
-(B::Blade{<:AbstractFloat}) = Blade(B, volume=-volume(B))
-(B::Pseudoscalar) = Pseudoscalar(B, value=-value(B))

"""
    reciprocal(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Scalar) = Scalar(B, value=1 / value(B))

reciprocal(B::Blade) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / norm(B)) :
        Blade(B, volume=-1 / norm(B))

reciprocal(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Pseudoscalar(B, value=1 / value(B)) :
        Pseudoscalar(B, value=-1 / value(B))

"""
    reverse(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reverse(B::Scalar) = B

reverse(B::Blade) =
    mod(grade(B), 4) < 2 ?  B : Blade(B, value=-value(B))

reverse(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?  B : Pseudoscalar(B, value=-value(B))


# --- Binary operations

# Imports
import Base.:(*)
import LinearAlgebra.:(⋅), LinearAlgebra.dot

# Exports
export ⋅, dot

"""
    ⋅(B::Union{AbstractBlade, Vector, Real},
      C::Union{AbstractBlade, Vector, Real})

    inner(B::Union{AbstractBlade, Vector, Real},
          C::Union{AbstractBlade, Vector, Real})

Return the inner product (left contraction) of `B` and `C`.
"""
# Inner products involving Scalars
⋅(B::Scalar, C::Scalar) = B * C
⋅(B::Scalar, C::Blade) = B * C
⋅(B::Blade, C::Scalar) = B * C
⋅(x::Real, B::Blade) = x * B
⋅(B::Blade, x::Real) = B * x
⋅(B::Scalar, C::Pseudoscalar) = B * C
⋅(B::Pseudoscalar, C::Scalar) = B * C

# Inner products involving Pseudoscalars
⋅(B::Pseudoscalar, C::Blade) = zero(B)
⋅(B::Blade, C::Pseudoscalar) = dual(B, C)

⋅(B::Pseudoscalar, C::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))

# Inner product between Blades
function ⋅(B::Blade{<:AbstractFloat}, C::Blade{<:AbstractFloat})
    assert_dim_equal(B, C)

    # (B ⋅ C) = 0 if grade(B) > grade(C)
    if grade(B) > grade(C)
        return zero(B)
    end

    # --- Compute (B ⋅ C) = proj(B, C) * C = proj(B, C) / C

    # Compute proj(B, C) = (B ⋅ C) / C
    projection = project(B, C)

    # Compute volume(B ⋅ C) = ±(norm(proj(B, C) * volume(C))
    volume_B_dot_C = (mod(grade(C), 4) < 2) ?
        volume(projection) * volume(C) :
       -volume(projection) * volume(C)

   # Construct blade representing dual of proj(B, C) scaled by volume(C)
    Blade(dual(projection, C), volume=volume_B_dot_C, copy_basis=false)
end

# Inner products between vectors and Blades
function ⋅(v::Vector{<:Real}, B::Blade)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    # --- Compute (v ⋅ B) = proj(v, B) * B = proj(v, B) / B

    # Compute proj(v, B) = (v ⋅ B) / B
    projection = project(v, B, return_blade=false)

    # Compute volume(v ⋅ B) = ±(norm(proj(v, B) * volume(B))
    volume_v_dot_B = (mod(grade(B), 4) < 2) ?
        LinearAlgebra.norm(projection) * volume(B) :
       -LinearAlgebra.norm(projection) * volume(B)

   # Construct blade representing dual of proj(v, B) scaled by volume(B)
    Blade(dual(Blade(projection), B), volume=volume_v_dot_B, copy_basis=false)
end

⋅(B::Blade, v::Vector{<:Real}) =
    length(v) != dim(B) ?
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`")) :
        dim(B) > 1 ? zero(B) : basis(B) ⋅ v

# Function aliases
dot(x::AbstractBlade, y::AbstractBlade) = x ⋅ y
dot(x::AbstractBlade, y::Union{Vector{<:Real}, Real}) = x ⋅ y
dot(x::Union{Vector{<:Real}, Real}, y::AbstractBlade) = x ⋅ y

"""
    *(B::Union{AbstractBlade}, C::Union{AbstractBlade})

Return the geometric product of `B` and `C`.
"""
# Geometric products involving Pseudoscalars
*(B::Pseudoscalar, C::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))

*(B::Blade, C::Pseudoscalar) = B ⋅ C
*(B::Pseudoscalar, C::Blade) = zero(B)

# Geometric product between Blades
function *(B::Blade{<:AbstractFloat}, C::Blade{<:AbstractFloat})
    assert_dim_equal(B, C)
    nothing  # TODO: implement
end

# Geometric product between vectors Blades
function *(v::Vector{<:Real}, B::Blade)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end
    v_dot_B = v ⋅ B
    v_wedge_B = v ∧ B

    if v_dot_B == zero(B)
        return v_wedge_B
    elseif v_wedge_B == zero(B)
        return v_dot_B
    end

    Multivector([v_dot_B, v_wedge_B])
end

*(B::Blade, v::Vector{<:Real}) = B ∧ v

*(v::Vector{<:Real}, w::Vector{<:Real}) =
    length(v) != length(w) ?
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`")) :
        v ⋅ w


# --- Utility functions

export rejection

"""
    rejection(vectors::Matrix, B::Blade; normalize::Bool=false)

Compute rejections of `vectors` from `B`. When `normalize` is true, the
rejection vectors are normalized.
"""
function rejection(vectors::Matrix, B::Blade; normalize::Bool=false)
    # --- Validate arguments

    if size(vectors, 1) != dim(B)
        throw(DimensionMismatch("`dim(vectors)` not equal to `dim(B)`"))
    end

    # --- Compute rejections using modified Gram-Schmidt algorithm

    # Initialize rejections
    rejections = Matrix{typeof(volume(B))}(vectors)

    # Remove basis(B) from rejections
    for idx_B in 1:grade(B)
        B_column = basis(B)[:, idx_B]
        for idx in 1:size(rejections, 2)
            rejections[:, idx] -=
                (rejections[:, idx] ⋅ B_column) * B_column
        end
    end

    # Normalize rejection vectors
    if normalize
        for idx in 1:size(rejections, 2)
            norm = LinearAlgebra.norm(rejections[:, idx])
            if norm > 0
                rejections[:, idx] /= norm
            end
        end
    end

    return rejections
end


# --- Non-exported utility functions

"""
    dim_equal(B::AbstractBlade, C::AbstractBlade)

Assert that the dimensions of `B` and `C` are the same.
"""
function assert_dim_equal(B::AbstractBlade, C::AbstractBlade)
    if dim(B) != dim(C)
        throw(DimensionMismatch("`dim(B)` not equal to `dim(C)`"))
    end
end
