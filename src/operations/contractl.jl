"""
contractl.jl defines methods for the contractl(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

export contractl
import Base.:(<)

# --- Operator aliases

<(M::AbstractMultivector, N::AbstractMultivector) = contractl(M, N)

<(M::AbstractMultivector, x::Real) = contractl(M, x)
<(x::Real, M::AbstractMultivector) = contractl(x, M)

<(M::AbstractMultivector, v::Vector{<:Real}) = contractl(M, v)
<(v::Vector{<:Real}, M::AbstractMultivector) = contractl(v, M)

# --- Method definitions

"""
    contractl(M, N)
    M < N

Compute the left contraction of the multivector `M` with the multivector `N`.
"""
contractl(M::AbstractMultivector, N::AbstractMultivector) = nothing  # TODO

# --- Operations involving an AbstractMultivector instance

# M::AbstractMultivector, B::AbstractBlade
# B::AbstractBlade, M::AbstractMultivector
contractl(M::AbstractMultivector, B::AbstractBlade) =
    Multivector(map(C -> contractl(C, B), blades(M)))

contractl(B::AbstractBlade, M::AbstractMultivector) =
    Multivector(map(C -> contractl(B, C), blades(M)))

# M::AbstractMultivector, B::AbstractScalar
# B::AbstractScalar, M::AbstractMultivector
contractl(M::AbstractMultivector, B::AbstractScalar) =
    length(M[0]) > 0 ?
        Scalar(value(M[0][1]) * B) :
        zero(Scalar{typeof(norm(M))})

contractl(B::AbstractScalar, M::AbstractMultivector) = B * M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
contractl(B::Zero, M::AbstractMultivector) = B
contractl(M::AbstractMultivector, B::Zero) = B

# M::AbstractMultivector, x::Real
# x::Real, M::AbstractMultivector
contractl(M::AbstractMultivector, x::Real) =
    length(M[0]) > 0 ?
        Scalar(value(M[0][1]) * x) :
        zero(Scalar{typeof(norm(M))})

contractl(x::Real, M::AbstractMultivector) = x * M

# --- Operations involving an AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
contractl(B::AbstractBlade, C::AbstractScalar) = zero(B)

# B::AbstractScalar, C::AbstractBlade
contractl(B::AbstractScalar, C::AbstractBlade) =
    Blade(C, volume=value(B) * volume(C))

# --- Operations involving a Blade instance

# B::Blade, C::Blade
function contractl(B::Blade, C::Blade)
    # --- Check arguments

    assert_dim_equal(B, C)

    # --- Handle edge cases

    # (B ⋅ C) = 0 if grade(B) > grade(C)
    if grade(B) > grade(C)
        return zero(B)
    end

    # --- Compute (B ⋅ C) = project(B, C) * C
    #     = volume(C) dual(project(B, C), C) * I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) project(B, C) / I_C
    #     = (-1)^((grade(C) * grade(C) - 1) / 2) volume(C)
    #       dual(project(B, C), C)
    #
    #     where I_C is the unit blade for the subspace represented by blade `C`
    #     that has the same orientation as basis(C).

    # Compute project(B, C) = (B ⋅ C) / C
    projection = project(B, C)

    # Compute
    #   (-1)^((grade(C) * grade(C) - 1) / 2) volume(C) dual(project(B, C), C)
    mod(grade(C), 4) < 2 ?
        volume(C) * dual(projection, C) :
       -volume(C) * dual(projection, C)
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function contractl(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    mod(grade(C), 4) < 2 ?
        dual(B, C) * volume(C) :
       -dual(B, C) * volume(C)
end

function contractl(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

# B::Blade, v::Vector
# v::Vector, B::Blade
function contractl(B::Blade, v::Vector{<:Real})
    assert_dim_equal(v, B)
    grade(B) > 1 ? zero(B) : Scalar(volume(B) * basis(B) ⋅ v)
end

function contractl(v::Vector{<:Real}, B::Blade)
    # --- Check arguments

    assert_dim_equal(v, B)

    # --- Compute (v ⋅ B) = project(v, B) * B
    #     = volume(B) dual(project(v, B), B) * I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) project(v, B) / I_B
    #     = (-1)^((grade(B) * grade(B) - 1) / 2) volume(B)
    #       dual(project(v, B), B)
    #
    #     where I_B is the unit blade for the subspace represented by blade `B`
    #     that has the same orientation as basis(B).

    # Compute project(v, B) = (v ⋅ B) / B
    projection = project(v, B)

    # Compute
    #   (-1)^((grade(B) * grade(B) - 1) / 2) volume(B) dual(project(v, B), B)
    mod(grade(B), 4) < 2 ?
        volume(B) * dual(projection, B) :
       -volume(B) * dual(projection, B)
end

# B::Blade, x::Real
# x::Real, B::Blade
contractl(B::Blade, x::Real) = zero(B)
contractl(x::Real, B::Blade) = x * B

# --- Operations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function contractl(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    mod(grade(B), 4) < 2 ?
        Scalar(value(B) * value(C)) :
        Scalar(-value(B) * value(C))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, B::Pseudoscalar
contractl(B::Pseudoscalar, C::AbstractScalar) = zero(B)
contractl(B::AbstractScalar, C::Pseudoscalar) =
    Pseudoscalar(C, value=value(B) * value(C))

# B::Pseudoscalar, v::Vector
# v::Vector, B::Pseudoscalar
function contractl(B::Pseudoscalar, v::Vector{<:Real})
    assert_dim_equal(B, v)
    zero(B)
end

function contractl(v::Vector{<:Real}, B::Pseudoscalar)
    assert_dim_equal(B, v)

    mod(grade(B), 4) < 2 ?
        volume(B) * dual(Blade(v)) :
       -volume(B) * dual(Blade(v))
end

# --- Operations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
contractl(B::AbstractScalar, C::One) = B
contractl(B::One, C::AbstractScalar) = C

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
contractl(B::AbstractScalar, C::Zero) = C
contractl(B::Zero, C::AbstractScalar) = B

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
contractl(B::AbstractScalar, x::Real) = B * x
contractl(x::Real, B::AbstractScalar) = x * B

# B::AbstractScalar, v::Vector
# v::Vector, B::AbstractScalar
contractl(B::AbstractScalar, v::Vector{<:Real}) = Blade(value(B) * v)
contractl(v::Vector{<:Real}, B::AbstractScalar) = zero(B)

# --- Operations involving a Scalar instance

# B::Scalar, C::Scalar
contractl(B::Scalar, C::Scalar) = B * C

# --- Operations involving a One instance

# B::One, C::One
contractl(B::One, C::One) = B

# B::One, C::Zero
# B::Zero, C::One
contractl(B::One, C::Zero) = C
contractl(B::Zero, C::One) = B

# --- Operations involving a Zero instance

# B::Zero, C::Zero
contractl(B::Zero, C::Zero) = B
