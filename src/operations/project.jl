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
project.jl defines methods for the project(x, y) function
"""

# --- Exports

export project

# --- Method definitions

using LinearAlgebra: I

"""
    project(M::AbstractMultivector, B::AbstractBlade)::AbstractMultivector

Compute the projection of `M` onto the subspace represented by `B`.

    project(B::AbstractBlade, C::AbstractBlade;
            return_blade=true)::Union{AbstractBlade, AbstractFloat,
                                      Vector, Matrix, LinearAlgebra.I}

Compute the projection of `B` onto the subspace represented by `C`.

When `return_blade` is true, the return value is an AbstractBlade. Otherwise,
the return value is an AbstractFloat if the result is a scalar, a Vector if the
result is a vector, a Matrix if the result is a blade with 1 < grade < `dim`,
and a multiple of LinearAlgebra.I if the result is a pseudoscalar.
"""
function project end

# ------ Specializations involving an AbstractMultivector instance

# M::AbstractMultivector, B::Zero
project(M::AbstractMultivector, B::Zero) = B

# ------ Specializations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
project(B::AbstractScalar, C::AbstractBlade; return_blade::Bool=true) =
    return_blade ? B : value(B)
project(B::AbstractBlade, C::AbstractScalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# B::AbstractBlade, C::Zero
project(B::AbstractBlade, C::Zero; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# B::AbstractBlade, C::Real
# B::Real, C::AbstractBlade
project(B::AbstractBlade, C::Real; return_blade::Bool=true) =
    return_blade ? zero(B) : 0
project(B::Real, C::AbstractBlade; return_blade::Bool=true) =
    return_blade ? Scalar(B) : B

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function project(B::Blade, C::Blade; return_blade::Bool=true)
    # --- Check arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # When grade(B) > grade(C), the projection is zero.
    if grade(B) > grade(C)
        return return_blade ? zero(B) : 0
    end

    # Compute the projections of basis(B) onto basis(C)
    projections = Matrix{typeof(volume(B))}(undef, dim(B), grade(B))
    for i in 1:grade(B)
        projections[:, i] = project(basis(B)[:, i], C, return_blade=false)
    end

    # Incorporate volume(B) into the norm of the first projection vector
    projections[:, 1] *= volume(B)

    # Compute projection (using fact that projection is an outermorphism)
    return_blade ? Blade(projections) : projections
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function project(B::Blade, C::Pseudoscalar; return_blade::Bool=true)
    assert_dim_equal(B, C)
    if return_blade
        return B
    else
        projections = copy(basis(B))
        projections[:, 1] *= volume(B)
        return projections
    end
end

function project(B::Pseudoscalar, C::Blade; return_blade::Bool=true)
    assert_dim_equal(B, C)
    return_blade ? zero(B) : 0
end

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

    return_blade ? Blade(projection) : vec(projection)
end

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function project(B::Pseudoscalar, C::Pseudoscalar; return_blade=true)
    assert_dim_equal(B, C)
    return_blade ? B : value(B) * I
end

# B::Pseudoscalar, v::AbstractScalar
project(B::Pseudoscalar, C::AbstractScalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# B::Pseudoscalar, C::Zero
project(B::Pseudoscalar, C::Zero; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

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

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::AbstractScalar
project(B::AbstractScalar, C::AbstractScalar; return_blade::Bool=true) =
    return_blade ? B : value(B)

# B::AbstractScalar, C::Zero
project(B::AbstractScalar, C::Zero; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
project(B::AbstractScalar, x::Real; return_blade::Bool=true) =
    return_blade ? B : value(B)

project(x::Real, B::AbstractScalar; return_blade::Bool=true) =
    return_blade ? Scalar{typeof(value(B))}(x) : x

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
project(B::AbstractScalar, v::Vector{<:Real}; return_blade::Bool=true) =
    return_blade ? B : value(B)

project(v::Vector{<:Real}, B::AbstractScalar; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# ------ Specializations involving a Zero instance

# B::Zero, x::Real
# x::Real, B::Zero
project(B::Zero, x::Real; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

project(x::Real, B::Zero; return_blade::Bool=true) =
    return_blade ? zero(B) : 0

# B::Zero, v::Vector
# v::Vector, B::Zero
project(B::Zero, v::Vector{<:Real}) = B
project(v::Vector{<:Real}, B::Zero) = B
