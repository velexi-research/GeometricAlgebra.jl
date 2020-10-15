"""
The operators.jl submodule defines operations on subtypes of AbstractBlade.

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
*(x::Real, B::Scalar) = Scalar(x * value(B))
*(B::Scalar, x::Real) = x * B
*(B::Scalar, C::Scalar) = Scalar(value(B) * value(C))

*(x::Real, B::Blade) = Blade(B, volume=x * volume(B))
*(B::Blade, x::Real) = x * B
*(B::Scalar, C::Blade) = value(B) * C
*(B::Blade, C::Scalar) = C * B

*(x::Real, B::Pseudoscalar) = Pseudoscalar(dim(B), x * value(B))
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
    project(v::Vector, B::AbstractBlade)
    project(B::AbstractBlade, v::Vector)

Return the projection of vector `v` onto the subspace represented by blade `B`.

When `B` is a Blade or a Pseudoscalar, the return value is a Vector. When `B`
is a Scalar, the return value is a Scalar representing zero.
"""
project(v::Vector{<:Real}, B::Scalar) = zero(B)
project(B::Scalar, v::Vector{<:Real}) = B

function project(v::Vector{<:Real}, B::Pseudoscalar)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end
    v
end

function project(B::Pseudoscalar, v::Vector{<:Real})
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end
    zero(B)
end

function project(v::Vector{<:Real}, B::Blade)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    grade(B) == 1 ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) :
        basis(B) * transpose(transpose(v) * basis(B))
end

function project(B::Blade, v::Vector{<:Real})
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    grade(B) == 1 ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) : zero(B)
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

    # Compute projection (using fact that projection is an outermorphism)
    projections = Matrix{typeof(volume(B))}(undef, dim(B), grade(B))
    for i in 1:grade(B)
        projections[:, i] = project(basis(B)[:, i], C)
    end
    Blade(projections)
end

"""
    dual(B::AbstractBlade)

Return the dual `B` (relative to the space that the geometric algebra is
extended from).
"""
dual(B::Scalar) = error("The dual of a Scalar is not well-defined")
dual(B::Pseudoscalar) = Scalar(value(B))

function dual(B::Blade)
    F = LinearAlgebra.qr(basis(B))
    dual_sign = mod(grade(B), 4) < 2 ?
        sign(LinearAlgebra.det(F.Q)) :
        -sign(LinearAlgebra.det(F.Q))

    dual_volume = mod(dim(B), 4) < 2 ?
        dual_sign * volume(B) :
        -dual_sign * volume(B)

    Blade(F.Q[:, grade(B) + 1:end], volume=dual_volume)
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
dual(B::Scalar, C::Blade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=value(B)) :
        Blade(C, volume=-value(B))

dual(B::Scalar, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=value(B)) :
        Pseudoscalar(C, value=-value(B))

# Duals involving Pseudoscalars
function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar(value(B))
end

function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

dual(B::Blade, C::Pseudoscalar) = dual(B)

# Duals involving Blades
function dual(B::Blade, C::Blade)
    # --- Validate arguments

    # Check that B and C are extended from the real vector spaces of the same
    # dimension
    assert_dim_equal(B, C)

    # Check that B is contained in C
    projection_coefficients = transpose(basis(B)) * basis(C)
    if LinearAlgebra.norm(projection_coefficients)^2 ≉ grade(B)
        error("`B` not contained in `C`")
    end

    # --- Handle edge cases

    # Subspaces represented by B and C are the same
    if grade(B) == grade(C)
        dual_sign = mod(grade(B), 4) < 2 ? 1 : -1

        dual_volume = LinearAlgebra.det(projection_coefficients) > 0 ?
            dual_sign * volume(B) : -dual_sign * volume(B)

        return Scalar(dual_volume)
    end

    # --- Compute dual using the basis vectors of `C` with the smallest
    #     projections (i.e., largest rejections) onto the basis of `B`.

    permutation = sortperm(sum(abs.(projection_coefficients), dims=1)[1, :])
    B_ext = hcat(basis(B), basis(C)[:, permutation[1:grade(C) - grade(B)]])
    F = LinearAlgebra.qr(B_ext)

    relative_sign_B_C = sign(prod(LinearAlgebra.diag(F.R)))
    dual_volume = mod(grade(B), 4) < 2 ?
        relative_sign_B_C * volume(B) :
        -relative_sign_B_C * volume(B)

    Blade(Matrix(F.Q)[:, grade(B) + 1:end], volume=dual_volume)
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
-(B::Scalar) = Scalar(-value(B))
-(B::Blade{<:AbstractFloat}) = Blade(B, volume=-volume(B))
-(B::Pseudoscalar) = Pseudoscalar(dim(B), -value(B))

"""
    reciprocal(B::AbstractBlade)

Return the multiplicative inverse of `B`.
"""
reciprocal(B::Scalar) = Scalar(1 / value(B))

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
# TODO: implement


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
⋅(B::Pseudoscalar, C::Pseudoscalar) = B * C
⋅(B::Pseudoscalar, C::Blade) = zero(B)
⋅(B::Blade, C::Pseudoscalar) = dual(B, C)

# Inner product between Blades
function ⋅(B::Blade{<:AbstractFloat}, C::Blade{<:AbstractFloat})
    assert_dim_equal(B, C)

    if grade(B) > grade(C)
        return zero(B)
    end

    nothing  # TODO: implement
end

# Inner products between vectors and Blades
function ⋅(v::Vector{<:Real}, B::Blade)
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    basis_B = basis(B)
    projection = zero(v)
    projection_norm = 0
    for k in 1:grade(B)
        projection_coefficient = v ⋅ basis_B[:, k]
        projection_norm += projection_coefficient ^ 2
        projection = projection_coefficient * basis_B[:, k]
    end
    projection_norm = sqrt(projection_norm)

    # TODO: fix sign
    Blade(dual(Blade(projection), B), volume=projection_norm * volume(B))
end

function ⋅(B::Blade, v::Vector{<:Real})
    if length(v) != dim(B)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    if dim(B) > 1
        return zero(B)
    end

    basis(B) ⋅ v
end

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

function *(v::Vector{<:Real}, w::Vector{<:Real})
    if length(v) != length(w)
        throw(DimensionMismatch("`dim(v)` not equal to `dim(B)`"))
    end

    v ⋅ w
end


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
