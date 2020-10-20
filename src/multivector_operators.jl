"""
The multivectcor_operators.jl submodule defines operations on subtypes of
AbstractMultivector.

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


# --- Core Multivector operations

export ∧, outer
export project, dual

# Note: scalar multiplication is grouped with geometric product functions.

"""
    ∧(M::AbstractMultivector, N::AbstractMultivector)
    ∧(M::AbstractMultivector, B::AbstractBlade)
    ∧(B::AbstractBlade, M::AbstractMultivector)
    ∧(M::AbstractMultivector, v::Vector)
    ∧(v::Vector, M::AbstractMultivector)
    ∧(M::AbstractMultivector, x::Real)
    ∧(x::Real, M::AbstractMultivector)

    outer(TODO)

Return the outer product of the arguments.
"""
# TODO


# --- Binary operations

# Imports
import Base.:(*)
import LinearAlgebra.:(⋅), LinearAlgebra.dot

"""
    ⋅(M::AbstractMultivector, N::AbstractMultivector)
    ⋅(M::AbstractMultivector, B::AbstractBlade)
    ⋅(B::AbstractBlade, M::AbstractMultivector)
    ⋅(M::AbstractMultivector, v::Vector)
    ⋅(v::Vector, M::AbstractMultivector)
    ⋅(M::AbstractMultivector, x::Real)
    ⋅(x::Real, M::AbstractMultivector)

    dot(TODO)

Return the inner product (left contraction) of the first argument with the
second argument.
"""
# TODO

"""
    *(M::AbstractMultivector, N::AbstractMultivector)
    *(M::AbstractMultivector, B::AbstractBlade)
    *(B::AbstractBlade, M::AbstractMultivector)
    *(M::AbstractMultivector, v::Vector)
    *(v::Vector, M::AbstractMultivector)

Return the geometric product of the arguments.

    *(M::AbstractMultivector, x::Real)
    *(x::Real, M::AbstractMultivector)

Return the product of the arguments when one (or both) of the arguments is a
scalar (i.e., scalar multiplication).
"""
# Geometric products involving scalars (i.e., scalar multiplication)
function *(x::Real, M::AbstractMultivector)
    blades = Vector{AbstractBlade}()
    for grade in grades(M)
        for B in M[grade]
        # for B in getindex(M, grade)
            push!(blades, x * B)
        end
    end
    Multivector(blades)
end
*(M::Multivector, x::Real) = x * M
*(B::Scalar, M::Multivector) = value(B) * M
*(M::Multivector, B::Scalar) = B * M
