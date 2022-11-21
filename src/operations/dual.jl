#   Copyright 2020 Velexi Corporation
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
dual.jl defines methods for the dual(x, y) function
"""

# --- Method definitions

# ------ Specializations involving a AbstractBlade instance

# B::AbstractBlade, C::AbstractScalar
# B::AbstractScalar, C::AbstractBlade
function dual(B::AbstractBlade, C::AbstractScalar)
    return grade(B) > 0 ? throw(ArgumentError("`B` not contained in `C`")) : B
end

function dual(B::AbstractScalar, C::AbstractBlade)
    return if mod(grade(C), 4) < 2
        Blade(C; volume=value(B), copy_basis=false)
    else
        Blade(C; volume=-value(B), copy_basis=false)
    end
end

# B::AbstractBlade, C::Real
# B::Real, C::AbstractBlade
function dual(B::AbstractBlade, C::Real)
    return grade(B) > 0 ? throw(ArgumentError("`B` not contained in `C`")) : B
end

function dual(B::Real, C::AbstractBlade)
    return if mod(grade(C), 4) < 2
        Blade(C; volume=B, copy_basis=false)
    else
        Blade(C; volume=-B, copy_basis=false)
    end
end

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

        dual_volume = if det(projection_coefficients) > 0
            dual_sign * volume(B)
        else
            -dual_sign * volume(B)
        end

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

    return Blade{typeof(volume(B))}(
        dim(B),
        grade(C) - grade(B),
        basis(C) * F.Q[:, (grade(B) + 1):end],
        dual_volume;
        enforce_constraints=false,
        copy_basis=false,
    )
end

# B::Blade, C::Pseudoscalar
# B::Pseudoscalar, C::Blade
function dual(B::Blade, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return dual(B)
end

function dual(B::Pseudoscalar, C::Blade)
    assert_dim_equal(B, C)
    throw(ArgumentError("`B` not contained in `C`"))
end

# B::Blade, C::Vector
# B::Vector, C::Blade
@inline dual(B::Blade, C::Vector) = dual(B, Blade(C))
@inline dual(B::Vector, C::Blade) = dual(Blade(B), C)

# ------ Specializations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
function dual(B::Pseudoscalar, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return Scalar{typeof(value(B))}(value(B))
end

# B::Pseudoscalar, C::AbstractScalar
# B::AbstractScalar, C::Pseudoscalar
dual(B::Pseudoscalar, C::AbstractScalar) = throw(ArgumentError("`B` not contained in `C`"))

function dual(B::AbstractScalar, C::Pseudoscalar)
    return if mod(grade(C), 4) < 2
        Pseudoscalar(C; value=value(B))
    else
        Pseudoscalar(C; value=-value(B))
    end
end

# B::Pseudoscalar, C::Real
# B::Real, C::Pseudoscalar
dual(B::Pseudoscalar, C::Real) = throw(ArgumentError("`B` not contained in `C`"))

function dual(B::Real, C::Pseudoscalar)
    return mod(grade(C), 4) < 2 ? Pseudoscalar(C; value=B) : Pseudoscalar(C; value=-B)
end

# B::Vector, C::Pseudoscalar
# B::Pseudoscalar, C::Vector
function dual(B::Vector, C::Pseudoscalar)
    assert_dim_equal(B, C)
    return dual(Blade(B))
end

function dual(B::Pseudoscalar, C::Vector)
    assert_dim_equal(B, C)

    if dim(B) == 1
        return dual(B, Blade(C))
    end

    throw(ArgumentError("`B` not contained in `C`"))
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

    return B
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

    return B
end

function dual(B::Real, C::AbstractScalar)
    if iszero(C)
        dual_relative_to_zero()
    end

    if iszero(B)
        dual_of_zero()
    end

    return Scalar(B)
end

# B::Vector, C::AbstractScalar
# B::AbstractScalar, C::Vector
function dual(B::Vector, C::AbstractScalar)
    return iszero(B) ? dual_of_zero() : throw(ArgumentError("`B` not contained in `C`"))
end
dual(B::AbstractScalar, C::Vector) = Blade(C; volume=value(B))

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
@inline function dual_relative_to_zero()
    return error("The dual of anything relative to Zero is not well-defined")
end

dual(M::Multivector, C::Zero) = dual_relative_to_zero()
dual(B::Blade, C::Zero) = dual_relative_to_zero()
dual(B::Pseudoscalar, C::Zero) = dual_relative_to_zero()
dual(B::Scalar, C::Zero) = dual_relative_to_zero()
dual(B::One, C::Zero) = dual_relative_to_zero()
dual(B::Zero, C::Zero) = dual_relative_to_zero()
dual(B::Real, C::Zero) = dual_relative_to_zero()
dual(B::Vector, C::Zero) = dual_relative_to_zero()
