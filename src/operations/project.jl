"""
project.jl defines methods for the project(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

export project

# --- Method definitions

"""
    project(M, B)

Compute the projection of multivector `M` onto the subspace represented by
blade `B`.
"""
project(M::AbstractMultivector, B::AbstractBlade) = nothing

# --- Operations involving an AbstractMultivector instance

# M::AbstractMultivector, B::Zero
project(M::AbstractMultivector, B::Zero) = B

# --- Operations involving an AbstractBlade instance

"""
    project(B, C; return_blade=true)

Compute the projection of blade `B` onto the subspace represented by blade `C`.

When `return_blade` is true, the return value is an AbstractBlade. Otherwise,
the return value is a Real (if the result is a scalar) or a Vector (if the
result is a vector).
"""
project(B::AbstractBlade, C::AbstractBlade; return_blade::Bool=true) = nothing

# B::AbstractScalar, C::AbstractBlade
project(B::AbstractScalar, C::AbstractBlade; return_blade::Bool=true) =
    return_blade ? B : value(B)

# B::AbstractBlade, v::Vector
# v::Vector, B::AbstractBlade
function project(B::AbstractBlade, v::Vector{<:Real}; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        v * LinearAlgebra.dot(v, basis(B)) : 0

    return_blade ? Blade(projection) : projection
end

function project(v::Vector{<:Real}, B::AbstractBlade; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) :
        basis(B) * transpose(transpose(v) * basis(B))

    return_blade ? Blade(projection) : projection
end

# --- Operations involving a Blade instance

# B::Blade, C::Blade
function project(B::Blade, C::Blade)
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
        projections[:, i] = project(basis(B)[:, i], C, return_blade=false)
    end

    # Incorporate volume(B) into the norm of the first projection vector
    projections[:, 1] *= volume(B)

    # Compute projection (using fact that projection is an outermorphism)
    Blade(projections)
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function project(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

function project(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

# B::Blade, C::AbstractScalar
project(B::Blade, C::AbstractScalar; return_blade::Bool=true) = zero(B)

# B::Blade, v::Vector
# v::Vector, B::Blade
function project(B::Blade, v::Vector{<:Real}; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        v * LinearAlgebra.dot(v, basis(B)) : 0

    return_blade ? Blade(projection) : projection
end

function project(v::Vector{<:Real}, B::Blade; return_blade::Bool=true)
    # Check arguments
    assert_dim_equal(v, B)

    # Compute projection
    projection = (grade(B) == 1) ?
        basis(B) * LinearAlgebra.dot(v, basis(B)) :
        basis(B) * transpose(transpose(v) * basis(B))

    return_blade ? Blade(projection) : projection
end

# --- Operations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function project(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    B
end

# B::Pseudoscalar, v::AbstractScalar
project(B::Pseudoscalar, C::AbstractScalar; return_blade::Bool=true) = zero(B)

# B::Pseudoscalar, v::Vector
# v::Vector, B::Pseudoscalar
function project(B::Pseudoscalar, v::Vector{<:Real}; return_blade::Bool=true)
    assert_dim_equal(v, B)
    return_blade ? zero(B) : 0
end

function project(v::Vector{<:Real}, B::Pseudoscalar; return_blade::Bool=true)
    assert_dim_equal(v, B)
    return_blade ? Blade(v) : v
end

# --- Operations involving an AbstractScalar instance

# B::AbstractScalar, C::AbstractScalar
project(B::AbstractScalar, C::AbstractScalar) = B

# B::AbstractScalar, C::Zero
project(B::AbstractScalar, C::Zero) = C

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
project(B::AbstractScalar, x::Real) = B
project(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x)

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
project(B::AbstractScalar, v::Vector{<:Real}; return_blade::Bool=true) =
    return_blade ? B : value(B)

project(v::Vector{<:Real}, B::AbstractScalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# --- Operations involving a Zero instance

# B::Zero, x::Real
# x::Real, B::Zero
project(B::Zero, x::Real) = B
project(x::Real, B::Zero) = B

# B::Zero, v::Vector
# v::Vector, B::Zero
project(B::Zero, v::Vector{<:Real}) = B
project(v::Vector{<:Real}, B::Zero) = B
