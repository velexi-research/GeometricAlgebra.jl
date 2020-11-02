"""
add.jl defines methods for the +(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

import Base.:(+)

# --- Method definitions

"""
    M + N

Compute the sum of the multivectors `M` and `N`.
"""
+(M::AbstractMultivector, N::AbstractMultivector) = nothing

# --- Operations involving an AbstractMultivector instance

# M::AbstractMultivector, B::One
# B::One,M::AbstractMultivector
+(M::AbstractMultivector, B::One) = Multivector(vcat([B], blades(M)))
+(B::One, M::AbstractMultivector) = B + M

# M::AbstractMultivector, B::Zero
# B::Zero, M::AbstractMultivector
+(M::AbstractMultivector, B::Zero) = M
+(B::Zero, M::AbstractMultivector) = M

# --- Operations involving a Blade instance

# B::Blade, C::Blade
function +(B::Blade, C::Blade)
    assert_dim_equal(B, C)

    if grade(B) == grade(C) == 1
        return Blade(volume(B) * basis(B) + volume(C) * basis(C))
    end

    # TODO: general blade sum
end

# --- Operations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
+(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) + value(C))

# --- Operations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
+(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) + 1)
+(B::One, C::AbstractScalar) = C + B

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
+(B::AbstractScalar, C::Zero) = B
+(B::Zero, C::AbstractScalar) = C

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
+(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) + x)
+(x::Real, B::AbstractScalar) = B + x

# --- Operations involving a Scalar instance

# B::Scalar, C::Scalar
+(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) + value(C))

# --- Operations involving a One instance

# B::One, C::One
+(B::One, C::One) = Scalar{typeof(value(B))}(2)

# B::One, C::Zero
# B::Zero, C::One
+(B::One, C::Zero) = B
+(B::Zero, C::One) = C

# --- Operations involving a Zero instance

# B::Zero, C::Zero
+(B::Zero, C::Zero) = B