"""
AbstractMultivector_operators.jl defines operators for the AbstractMultivector
type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Comparison operators

# Note: the comparison operators defined here ensure that arbitrary
#       AbstractMultivector instances can be compared

"""
    ==(M::AbstractMultivector, N::AbstractMultivector)

Return true if M and N are equal; otherwise, return false.
"""
==(B::AbstractMultivector, C::AbstractMultivector) = false

"""
    ≈(M::AbstractMultivector, N::AbstractMultivector)

Return true if `M` and `N` are approximately equal; otherwise, return false.
"""
≈(M::AbstractMultivector, N::AbstractMultivector) = false
≈(M::AbstractMultivector, x::Real) = false
≈(x::Real, M::AbstractMultivector) = false

# --- Unary operators

# TODO

# --- Binary operators

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

TODO: Review and revise
"""
# TODO

"""
    wedge(M, N)
    M ∧ N

Valid arguments
---------------
    ∧(M::AbstractMultivector, N::AbstractMultivector)
    ∧(M::AbstractMultivector, B::AbstractBlade)
    ∧(B::AbstractBlade, M::AbstractMultivector)
    ∧(M::AbstractMultivector, v::Vector)
    ∧(v::Vector, M::AbstractMultivector)
    ∧(M::AbstractMultivector, x::Real)
    ∧(x::Real, M::AbstractMultivector)

Return the outer product of the arguments.

TODO: Review and revise
"""
wedge(M::AbstractMultivector, N::AbstractMultivector) = nothing  # TODO

"""
    contractl(M, N)

TODO: Review and revise
"""
contractl(M::AbstractMultivector, N::AbstractMultivector) = nothing  # TODO

# --- Operator aliases

∧(M::AbstractMultivector, N::AbstractMultivector) = wedge(M, N)
# const ∧ = wedge

"""
    dot(M, N)
    M ⋅ N

Valid arguments
---------------
    dot(M::AbstractMultivector, N::AbstractMultivector)
    dot(M::AbstractMultivector, B::AbstractBlade)
    dot(B::AbstractBlade, M::AbstractMultivector)
    dot(M::AbstractMultivector, v::Vector)
    dot(v::Vector, M::AbstractMultivector)
    dot(M::AbstractMultivector, x::Real)
    dot(x::Real, M::AbstractMultivector)

Return the inner product (left contraction) of the first argument with the
second argument.

TODO: Review and revise
"""
dot(M::AbstractMultivector, N::AbstractMultivector; left=true) = contractl(M, N)

<(M::AbstractMultivector, N::AbstractMultivector) = contractl(M, N)

# --- Special cases

# Operations involving Vectors
*(M::AbstractMultivector, v::Vector) = Multivector(map(B -> B * v, blades(M)))
*(v::Vector, M::AbstractMultivector) = Multivector(map(B -> v * B, blades(M)))

# Operations involving AbstractBlades
*(M::AbstractMultivector, B::AbstractBlade) =
    Multivector(map(C -> C * B, blades(M)))

*(B::AbstractBlade, M::AbstractMultivector) =
    Multivector(map(C -> B * C, blades(M)))

# Operations involving Reals
*(x::Real, M::AbstractMultivector) = Multivector(map(B -> x * B, blades(M)))
*(M::AbstractMultivector, x::Real) = x * M

wedge(M::AbstractMultivector, x::Real) = x * M
wedge(x::Real, M::AbstractMultivector) = x * M

contractl(M::AbstractMultivector, x::Real) =
    length(M[0]) > 0 ?
        Scalar(value(M[0][1]) * x) :
        zero(Scalar{typeof(norm(M))})

contractl(x::Real, M::AbstractMultivector) = x * M

∧(M::AbstractMultivector, x::Real) = wedge(M, x)
∧(x::Real, M::AbstractMultivector) = wedge(x, M)

dot(M::AbstractMultivector, x::Real; left=true) = contractl(M, x)
dot(x::Real, M::AbstractMultivector; left=true) = contractl(x, M)

<(M::AbstractMultivector, x::Real) = contractl(M, x)
<(x::Real, M::AbstractMultivector) = contractl(x, M)
