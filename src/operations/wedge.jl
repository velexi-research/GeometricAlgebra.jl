"""
wedge.jl defines methods for the wedge(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

export wedge, ∧

# --- Method definitions

"""
    wedge(M, N)
    M ∧ N

Compute the outer product of the multivector `M` with the multivector `N`.
"""
wedge(M::AbstractMultivector, N::AbstractMultivector) = nothing  # TODO

# M::AbstractMultivector, B::AbstractBlade
# B::AbstractBlade, M::AbstractMultivector
wedge(M::AbstractMultivector, B::AbstractBlade) =
    Multivector(map(C -> wedge(C, B), blades(M)))

wedge(B::AbstractBlade, M::AbstractMultivector) =
    Multivector(map(C -> wedge(B, C), blades(M)))

# M::AbstractMultivector, B::AbstractScalar
# B::AbstractScalar, M::AbstractMultivector
wedge(B::AbstractMultivector, C::AbstractScalar) = B * C
wedge(B::AbstractScalar, C::AbstractMultivector) = B * C

# M::AbstractMultivector, B::One
# B::One, M::AbstractMultivector
wedge(M::AbstractMultivector, B::One) = M
wedge(B::One, M::AbstractMultivector) = M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
wedge(B::Zero, M::AbstractMultivector) = B
wedge(M::AbstractMultivector, B::Zero) = B

# M::AbstractMultivector, x::Real
# x::Real, M::AbstractMultivector
wedge(M::AbstractMultivector, x::Real) = x * M
wedge(x::Real, M::AbstractMultivector) = x * M

# --- Operations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
wedge(B::AbstractBlade, C::AbstractScalar) = B * C
wedge(B::AbstractScalar, C::AbstractBlade) = B * C

# B::AbstractBlade, x::Real
# x::Real, B::AbstractBlade
wedge(B::AbstractBlade, x::Real) = B * x
wedge(x::Real, B::AbstractBlade) = x * B

# --- Operator aliases

∧(M::AbstractMultivector, N::AbstractMultivector) = wedge(M, N)

∧(M::AbstractMultivector, x::Real) = wedge(M, x)
∧(x::Real, M::AbstractMultivector) = wedge(x, M)

∧(M::AbstractMultivector, v::Vector{<:Real}) = wedge(M, Blade(v))
∧(v::Vector{<:Real}, M::AbstractMultivector) = wedge(Blade(v), M)

# --- Operations involving a Blade instance

# B::Blade, C::Blade
function wedge(B::Blade, C::Blade)
    # --- Check arguments

    assert_dim_equal(B, C)

    if grade(B) + grade(C) > dim(B)
        return zero(B)
    end

    # --- Construct new blade

#    # Incorporate volume(B) and volume(C) into the norm of the first vector of
#    # basis(B)
#    basis_B = copy(basis(B))
#    basis_B[:, 1] *= volume(B) * volume(C)

    volume(B) * volume(C) * Blade(hcat(basis(B), basis(C)))
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function wedge(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    zero(B)
end

wedge(B::Pseudoscalar, C::Blade) = wedge(C, B)

# B::Blade, v::Vector
# v::Vector, B::Blade
function wedge(B::Blade, v::Vector{<:Real})
    assert_dim_equal(B, v)

    # Note: volume(B) is incorporated into the norm of `v`
    Blade(hcat(basis(B), volume(B) * v))
end

wedge(v::Vector{<:Real}, B::Blade) =
    mod(grade(B), 2) == 0 ?
        wedge(B, v) : wedge(B, -v)

# v::Vector, w::Vector
wedge(v::Vector{<:Real}, w::Vector{<:Real}) = Blade(hcat(v, w))

# --- Operations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function wedge(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    zero(B)
end

# B::Pseudoscalar, v::Vector
# v::Vector, B::Pseudoscalar
function wedge(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    zero(B)
end

wedge(v::Vector{<:Real}, B::Pseudoscalar) = B ∧ v

# --- Operations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
wedge(B::AbstractScalar, C::One) = B
wedge(B::One, C::AbstractScalar) = C

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
wedge(B::AbstractScalar, C::Zero) = C
wedge(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
wedge(B::AbstractScalar, x::Real) = B * x
wedge(x::Real, B::AbstractScalar) = x * B

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
wedge(B::AbstractScalar, v::Vector{<:Real}) =  Blade(value(B) * v)
wedge(v::Vector{<:Real}, B::AbstractScalar) = wedge(B, v)

# --- Operations involving a Scalar instance

# B::Scalar, C::Scalar
wedge(B::Scalar, C::Scalar) = B * C

# --- Operations involving a One instance

# B::One, C::One
wedge(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
wedge(B::One, C::Zero) = C
wedge(B::Zero, C::One) = B

# --- Operations involving a Zero instance

# B::Zero, C::Zero
wedge(B::Zero, C::Zero) = B
