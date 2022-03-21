"""
dual.jl defines methods for the dual(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Method definitions

# ------ Specializations involving a AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
dual(B::AbstractBlade, C::AbstractScalar) = grade(B) > 0 ? zero(B) : B

dual(B::AbstractScalar, C::AbstractBlade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=value(B), copy_basis=false) :
        Blade(C, volume=-value(B), copy_basis=false)

# B::AbstractBlade, C::Real
# B::Real, C::AbstractBlade
dual(B::AbstractBlade, C::Real) = grade(B) > 0 ? zero(B) : B

dual(B::Real, C::AbstractBlade) =
    mod(grade(C), 4) < 2 ?
        Blade(C, volume=B, copy_basis=false) :
        Blade(C, volume=-B, copy_basis=false)

# ------ Specializations involving a Blade instance

# B::Blade, C::Blade
function dual(B::Blade, C::Blade)
    # --- Check arguments

    # Check that B and C are from the same real vector space
    assert_dim_equal(B, C)

    # Check that B is contained in C
    if grade(B) > grade(C)
        throw(ArgumentError("`B` not contained in `C`"))
    end

    projection_coefficients = transpose(basis(C)) * basis(B)
    if LinearAlgebra.norm(projection_coefficients)^2 ≉ grade(B)
        throw(ArgumentError("`B` not contained in `C`"))
    end

    # --- Handle edge cases

    # Subspaces represented by B and C are the same
    if grade(B) == grade(C)
        dual_sign = mod(grade(B), 4) < 2 ? 1 : -1

        dual_volume = det(projection_coefficients) > 0 ?
            dual_sign * volume(B) : -dual_sign * volume(B)

        return Scalar{typeof(volume(B))}(dual_volume)
    end

    # --- Extend basis(B) to an orthonormal basis for entire subspace
    #     represented by basis(C)

    F = qr(projection_coefficients)

    # --- Compute volume of dual

    # Account for orientation of Q relative to orientation of basis(C)
    dual_volume = volume(B) * sign(det(F.Q))

    # Account for orientation of first grade(B) columns of Q relative to
    # orientation of basis(B)
    if prod(diag(F.R)) < 0
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

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function dual(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    dual(B)
end

function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    zero(B)
end

# B::Blade, C::Vector
# B::Vector, C::Blade
function dual(B::Blade, C::Vector)
    C_blade = Blade(C)
    return dual(B, C_blade)
end

function dual(B::Vector, C::Blade)
    B_blade = Blade(B)
    return dual(B_blade, C)
end

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    Scalar{typeof(value(B))}(value(B))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, C::Pseudoscalar
dual(B::Pseudoscalar, C::AbstractScalar) = zero(B)

dual(B::AbstractScalar, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=value(B)) :
        Pseudoscalar(C, value=-value(B))

# B::Pseudoscalar, C::Real
# B::Real, C::Pseudoscalar
dual(B::Pseudoscalar, C::Real) = zero(B)

dual(B::Real, C::Pseudoscalar) =
    mod(grade(C), 4) < 2 ?
        Pseudoscalar(C, value=B) :
        Pseudoscalar(C, value=-B)

# B::Vector, C::Pseudoscalar
# B::Pseudoscalar, C::Vector
function dual(B::Vector, C::Pseudoscalar)
    assert_dim_equal(B, C)
    dual(Blade(B))
end

function dual(B::Pseudoscalar, C::Vector)
    assert_dim_equal(B, C)

    if length(C) == 1
        return dual(B, Blade(C))
    end

    zero(B)
end

# ------ Specializations involving an AbstractScalar instance

# B::AbstractScalar, C::AbstractScalar
function dual(B::AbstractScalar, C::AbstractScalar)
    if iszero(C)
        dual_relative_to_zero()
    end

    if iszero(B)
        dual_of_zero()
    end

    B
end

# B::AbstractScalar, C::Real
# B::Real, C::AbstractScalar
function dual(B::AbstractScalar, C::Real)
    if iszero(C)
        dual_relative_to_zero()
    end

    if iszero(B)
        dual_of_zero()
    end

    B
end

function dual(B::Real, C::AbstractScalar)
    if iszero(C)
        dual_relative_to_zero()
    end

    if iszero(B)
        dual_of_zero()
    end

    Scalar(B)
end

# B::Vector, C::AbstractScalar
# B::AbstractScalar, C::Vector
dual(B::Vector, C::AbstractScalar) = zero(B)
dual(B::AbstractScalar, C::Vector) = Blade(C, volume=value(B))

# ------ Specializations involving a Zero instance

# dual(B::Zero, C)
@inline dual_of_zero() = error("The dual of Zero is not well-defined")

dual(B::Zero, C::Blade) = dual_of_zero()
dual(B::Zero, C::Pseudoscalar) = dual_of_zero()
dual(B::Zero, C::Scalar) = dual_of_zero()
dual(B::Zero, C::One) = dual_of_zero()
dual(B::Zero, C::Real) = dual_of_zero()
dual(B::Zero, C::Vector) = dual_of_zero()

# dual(B, C::Zero)
@inline dual_relative_to_zero() =
    error("The dual of anything relative to Zero is not well-defined")

dual(M::Multivector, C::Zero) = dual_relative_to_zero()
dual(B::Blade, C::Zero) = dual_relative_to_zero()
dual(B::Pseudoscalar, C::Zero) = dual_relative_to_zero()
dual(B::Scalar, C::Zero) = dual_relative_to_zero()
dual(B::One, C::Zero) = dual_relative_to_zero()
dual(B::Zero, C::Zero) = dual_relative_to_zero()
dual(B::Real, C::Zero) = dual_relative_to_zero()
dual(B::Vector, C::Zero) = dual_relative_to_zero()
