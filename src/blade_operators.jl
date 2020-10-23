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


# --- Core Blade operations

export wedge, ∧
export proj, dual

# Note: scalar multiplication is grouped with the geometric product functions.

"""
    wedge(B, C)
    B ∧ C

Return the outer product of the arguments.

Valid arguments
---------------
    wedge(B::AbstractBlade, C::AbstractBlade)
    wedge(B::AbstractBlade, v::Vector)
    wedge(v::Vector, B::AbstractBlade)
    wedge(v::Vector, w::Vector)
    wedge(B::AbstractBlade, x::Real)
    wedge(x::Real, B::AbstractBlade)
"""
wedge(B::Scalar, C::Scalar) = B * C
wedge(B::Scalar, C::Blade) = B * C
wedge(B::Blade, C::Scalar) = B * C
wedge(B::Scalar, C::Pseudoscalar) = B * C
wedge(B::Pseudoscalar, C::Scalar) = B * C

function wedge(B::Blade, C::Blade)
    assert_dim_equal(B, C)
    Blade(hcat(basis(B), basis(C)), volume=volume(B) * volume(C))
end

function wedge(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end
wedge(B::Blade, C::Pseudoscalar) = C ∧ B

function wedge(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    zero(B)
end

wedge(B::Scalar, v::Vector{<:Real}) =  Blade(value(B) * v)
wedge(v::Vector{<:Real}, B::Scalar) = B ∧ v

function wedge(v::Vector{<:Real}, B::Blade)
    assert_dim_equal(v, B)
    Blade(hcat(v, basis(B)), volume=LinearAlgebra.norm(v) * volume(B))
end

function wedge(B::Blade, v::Vector{<:Real})
    assert_dim_equal(B, v)
    Blade(hcat(basis(B), v), volume=LinearAlgebra.norm(v) * volume(B))
end

function wedge(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    zero(B)
end

wedge(v::Vector{<:Real}, B::Pseudoscalar) = B ∧ v

wedge(v::Vector{<:Real}, w::Vector{<:Real}) = Blade(hcat(v, w))

wedge(x::Real, B::Scalar) = x * B
wedge(B::Scalar, x::Real) = B * x
wedge(x::Real, B::Blade) = x * B
wedge(B::Blade, x::Real) = B * x
wedge(x::Real, B::Pseudoscalar) = x * B
wedge(B::Pseudoscalar, x::Real) = B * x

const ∧ = wedge

"""
    proj(B, C; return_blade=true)

Return the projection of blade `B` onto the subspace represented by blade `C`.

When `return_blade` is true, the return value is an AbstractBlade. Otherwise,
the return value is a Real (if the result is a scalar) or a Vector.

Valid arguments
---------------
    proj(B::AbstractBlade, C::AbstractBlade; return_blade::Bool)
    proj(v::Vector, B::AbstractBlade; return_blade::Bool)
    proj(B::AbstractBlade, v::Vector; return_blade::Bool)
"""
proj(v::Vector{<:Real}, B::Scalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

proj(B::Scalar, v::Vector{<:Real}; return_blade::Bool=true) =
    return_blade ? B : value(B)

function proj(v::Vector{<:Real}, B::Blade; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) :
        basis(B) * transpose(transpose(v) * basis(B))

    return_blade ? Blade(projection) : projection
end

function proj(B::Blade, v::Vector{<:Real}; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        v * LinearAlgebra.dot(v, basis(B)) : 0

    return_blade ? Blade(projection) : projection
end

function proj(v::Vector{<:Real}, B::Pseudoscalar; return_blade::Bool=true)
    assert_dim_equal(v, B)
    return_blade ? Blade(v) : v
end

function proj(B::Pseudoscalar, v::Vector{<:Real}; return_blade::Bool=true)
    assert_dim_equal(v, B)
    return_blade ? zero(B) : 0
end

proj(B::Scalar, C::Scalar) = B
proj(B::Scalar, C::Blade) = B
proj(B::Blade, C::Scalar) = zero(B)
proj(B::Scalar, C::Pseudoscalar) = B
proj(B::Pseudoscalar, C::Scalar) = zero(B)

function proj(B::Blade, C::Blade)
    # --- Check arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # When grade(B) > grade(C), the projection is zero.
    if grade(B) > grade(C)
        return zero(B)
    end

    # Compute the projections of basis(B) onto basis(C)
    projections = Matrix{typeof(volume(B))}(undef, dim(B), grade(B))
    for i in 1:grade(B)
        projections[:, i] = proj(basis(B)[:, i], C, return_blade=false)
    end

    # Encode volume(B) in the first projection vector
    projections[:, 1] *= volume(B)

    # Compute projection (using fact that projection is an outermorphism)
    Blade(projections)
end

function proj(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function proj(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

function proj(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

"""
    dual(B)

Return the dual `B` (relative to the space that the geometric algebra is
extended from).

Valid arguments
---------------
    dual(B::Blade)
    dual(B::Pseudoscalar)
"""
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

dual(B::Pseudoscalar) = Scalar(value(B))

"""
    dual(B, dim)

Return the dual of `B` when `B` is a Scalar. Note that the dimension of the
embedding space must be explicitly specified.

Valid arguments
---------------
    dual(B::Scalar, dim::Integer)
"""
dual(B::Scalar, dim::Integer) =
    mod(dim, 4) < 2 ?
        Pseudoscalar(dim, value(B)) :
        Pseudoscalar(dim, -value(B))

"""
    dual(B, C)

Return the dual `B` relative to the subspace represented by `C`.

Valid arguments
---------------
    dual(B::AbstractBlade, C::AbstractBlade)

Notes
-----
* `dual(B, C)` is only defined if (1) `B` and `C` are extended from real
  vector spaces of the same dimension and (2) the subspace represented by `B`
  is contained in subspace represented by `C`.

* The volume of `C` is ignored.
"""
dual(B::Scalar, C::Scalar) = B

dual(B::Scalar, C::Blade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=value(B), copy_basis=false) :
        Blade(C, volume=-value(B), copy_basis=false)

dual(B::Blade, C::Scalar) = zero(B)

dual(B::Scalar, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=value(B)) :
        Pseudoscalar(C, value=-value(B))

dual(B::Pseudoscalar, C::Scalar) = zero(B)

function dual(B::Blade, C::Blade)
    # --- Check arguments

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

    Blade{typeof(volume(B))}(dim(B), grade(C) - grade(B),
                             basis(C) * F.Q[:, grade(B) + 1:end], dual_volume,
                             enforce_constraints=false,
                             copy_basis=false)
end

function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function dual(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    dual(B)
end

function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar(value(B))
end


# --- Unary operations

# Imports
import Base.:(-), Base.reverse

# Exports
export reciprocal

"""
    -B

Return the additive inverse of `B`.

Valid arguments
---------------
    -(B::AbstractBlade)
"""
-(B::Scalar) = Scalar(B, value=-value(B))
-(B::Blade{<:AbstractFloat}) = Blade(B, volume=-volume(B), copy_basis=false)
-(B::Pseudoscalar) = Pseudoscalar(B, value=-value(B))

"""
    reciprocal(B)

Return the multiplicative inverse of `B`.

Valid arguments
---------------
    reciprocal(B::AbstractBlade)
"""
reciprocal(B::Scalar) = Scalar(B, value=1 / value(B))

reciprocal(B::Blade) =
    mod(grade(B), 4) < 2 ?
        Blade(B, volume=1 / volume(B), copy_basis=false) :
        Blade(B, volume=-1 / volume(B), copy_basis=false)

reciprocal(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?
        Pseudoscalar(B, value=1 / value(B)) :
        Pseudoscalar(B, value=-1 / value(B))

"""
    reverse(B)

Return the multiplicative inverse of `B`.

Valid arguments
---------------
    reverse(B::AbstractBlade)
"""
reverse(B::Scalar) = B

reverse(B::Blade) =
    mod(grade(B), 4) < 2 ?  B : Blade(B, volume=-volume(B), copy_basis=false)

reverse(B::Pseudoscalar) =
    mod(grade(B), 4) < 2 ?  B : Pseudoscalar(B, value=-value(B))


# --- Binary operations

# Imports
import Base.:(*)
import LinearAlgebra.dot

# Exports
export dot, ⋅

"""
    dot(B, C)
    B ⋅ C

Return the inner product (left contraction) of the first argument with the
second argument.

Valid arguments
---------------
    dot(B::AbstractBlade, C::AbstractBlade)
    dot(B::AbstractBlade, v::Vector)
    dot(v::Vector, B::AbstractBlade)
    dot(B::AbstractBlade, x::Real)
    dot(x::Real, B::AbstractBlade)
"""
dot(B::Scalar, C::Scalar) = B * C
dot(B::Scalar, C::Blade) = B * C
dot(B::Blade, C::Scalar) = zero(B)
dot(B::Scalar, C::Pseudoscalar) = B * C
dot(B::Pseudoscalar, C::Scalar) = zero(B)

function dot(B::Blade{<:Real}, C::Blade{<:Real})
    # --- Check arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # (B ⋅ C) = 0 if grade(B) > grade(C)
    if grade(B) > grade(C)
        return zero(B)
    end

    # --- Compute (B ⋅ C) = proj(B, C) * C
    #     = volume(C) dual(proj(B, C), C) * I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) proj(B, C) / I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) dual(proj(B, C), C)
    #
    #     where I_C is the unit blade for the subspace represented by blade `C`
    #     that has the same orientation as basis(C).

    # Compute proj(B, C) = (B ⋅ C) / C
    projection = proj(B, C)

    # Compute (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) dual(proj(B, C), C)
    mod(grade(C), 4) < 2 ?
        volume(C) * dual(projection, C) :
       -volume(C) * dual(projection, C)
end

function dot(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

function dot(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)

    mod(grade(C), 4) < 2 ?
        dual(B, C) * volume(C) :
       -dual(B, C) * volume(C)
end

function dot(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)

    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))
end

dot(B::Scalar, v::Vector{<:Real}) = B * v
dot(v::Vector{<:Real}, B::Scalar) = zero(B)

function dot(v::Vector{<:Real}, B::Blade)
    # --- Check arguments

    assert_dim_equal(v, B)

    # --- Compute (v ⋅ B) = proj(v, B) * B
    #     = volume(B) dual(proj(v, B), B) * I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) proj(v, B) / I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) dual(proj(v, B), B)
    #
    #     where I_B is the unit blade for the subspace represented by blade `B`
    #     that has the same orientation as basis(B).

    # Compute proj(v, B) = (v ⋅ B) / B
    projection = proj(v, B)

    # Compute (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) dual(proj(v, B), B)
    mod(grade(B), 4) < 2 ?
        volume(B) * dual(projection, B) :
       -volume(B) * dual(projection, B)
end

function dot(B::Blade, v::Vector{<:Real})
    assert_dim_equal(v, B)
    grade(B) > 1 ? zero(B) : Scalar(volume(B) * basis(B) ⋅ v)
end

function dot(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(v, B)
    zero(B)
end

function dot(v::Vector{<:Real}, B::Pseudoscalar)
    assert_dim_equal(v, B)
    mod(grade(B), 4) < 2 ?
        value(B) * dual(Blade(v)) :
       -value(B) * dual(Blade(v))
end

dot(x::Real, B::Scalar) = x * B
dot(B::Scalar, x::Real) = B * x
dot(x::Real, B::Blade) = x * B
dot(B::Blade, x::Real) = zero(B)
dot(x::Real, B::Pseudoscalar) = x * B
dot(B::Pseudoscalar, x::Real) = zero(B)

"""
    B * C

Return the geometric product of the arguments.
TODO: define a geometric product function

Valid arguments
---------------
    *(B::AbstractBlade, C::AbstractBlade)
    *(B::AbstractBlade, v::Vector)
    *(v::Vector, B::AbstractBlade)
    *(B::AbstractBlade, x::Real)
    *(x::Real, B::AbstractBlade)
"""
*(B::Scalar, C::Scalar) = Scalar(B, value=value(B) * value(C))
*(B::Scalar, C::Blade) = value(B) * C
*(B::Blade, C::Scalar) = C * B
*(B::Scalar, C::Pseudoscalar) = Pseudoscalar(C, value=value(B) * value(C))
*(B::Pseudoscalar, C::Scalar) = C * B

function *(B::Blade{<:AbstractFloat}, C::Blade{<:AbstractFloat})
    # --- Preparations

    # Check arguments
    assert_dim_equal(B, C)

    # --- Compute geometric product

    if grade(B) < grade(C)
        println("GOT HERE B")
        M = C
        for i in grade(B):-1:1
            M = basis(B)[:, i] * M
        end
        M = volume(B) * M
    else
        M = B
        for i in 1:grade(C)
            M = M * basis(C)[:, i]
        end
        M = volume(C) * M
    end

    M
end

*(B::Blade, C::Pseudoscalar) = B ⋅ C
*(B::Pseudoscalar, C::Blade) = zero(B)
*(B::Pseudoscalar, C::Pseudoscalar) = B ⋅ C

*(B::Scalar, v::Vector{<:Real}) = B * Blade(v)
*(v::Vector{<:Real}, B::Scalar) = Blade(v) * B
*(B::Pseudoscalar, v::Vector{<:Real}) = zero(B)
*(v::Vector{<:Real}, B::Pseudoscalar) = Blade(v) * B

function *(v::Vector{<:Real}, B::Blade)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute geometric product
    v_dot_B = v ⋅ B
    v_wedge_B = v ∧ B

    if v_dot_B == zero(v_dot_B)
        return v_wedge_B
    elseif v_wedge_B == zero(v_wedge_B)
        return v_dot_B
    end

    Multivector([v_dot_B, v_wedge_B])
end

function *(B::Blade, v::Vector{<:Real})
    # Check arguments
    assert_dim_equal(v, B)

    # Compute geometric product
    B_dot_v = B ⋅ v
    B_wedge_v = B ∧ v

    if B_dot_v == zero(B_dot_v)
        return B_wedge_v
    elseif B_wedge_v == zero(B_wedge_v)
        return B_dot_v
    end

    println(B_dot_v)
    println(B_wedge_v)
    Multivector([B_dot_v, B_wedge_v])
end

*(x::Real, B::Scalar) = Scalar(B, value=x * value(B))
*(B::Scalar, x::Real) = x * B
*(x::Real, B::Blade) = Blade(B, volume=x * volume(B), copy_basis=false)
*(B::Blade, x::Real) = x * B
*(x::Real, B::Pseudoscalar) = Pseudoscalar(B, value=x * value(B))
*(B::Pseudoscalar, x::Real) = x * B


# --- Utility functions

export rejection

"""
    rejection(vectors, B; normalize=false)

Compute rejections of `vectors` from `B`. When `normalize` is true, the
rejection vectors are normalized.

Valid arguments
---------------
    rejection(vectors::Matrix, B::Blade; normalize::Bool=false)
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
    dim_equal(B, C)

Assert that the dimensions of `B` and `C` are the same.

Valid arguments
---------------
    dim_equal(B::AbstractBlade, C::AbstractBlade)
    dim_equal(B::AbstractBlade, v::Vector)
    dim_equal(v::Vector, B::AbstractBlade)
"""
function assert_dim_equal(B::AbstractBlade, C::AbstractBlade)
    if dim(B) != dim(C)
        throw(DimensionMismatch("`dim(B)` not equal to `dim(C)`"))
    end
end

function assert_dim_equal(B::AbstractBlade, v::Vector)
    if dim(B) != length(v)
        throw(DimensionMismatch("`dim(B)` not equal to `length(v)`"))
    end
end
assert_dim_equal(v::Vector, B::AbstractBlade) = assert_dim_equal(B, v)
