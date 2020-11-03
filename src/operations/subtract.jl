"""
subtract.jl defines methods for the -(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

import Base.:(-)

# --- Method definitions

"""
    M - N

Compute the difference of the multivectors `M` and `N`.
"""
-(M::AbstractMultivector, N::AbstractMultivector) = M + -N

# --- Operations involving a Pseudoscalar instance

# B::Pseudoscalar, C::Pseudoscalar
-(B::Pseudoscalar, C::Pseudoscalar) = Pseudoscalar(B, value=value(B) - value(C))

# --- Operations involving an AbstractScalar instance

# B::AbstractScalar, C::One
# B::One, C::AbstractScalar
-(B::AbstractScalar, C::One) = Scalar{typeof(value(B))}(value(B) - 1)
-(B::One, C::AbstractScalar) = Scalar{typeof(value(C))}(1 - value(C))

# B::AbstractScalar, C::Zero
# B::Zero, C::AbstractScalar
-(B::AbstractScalar, C::Zero) = B
-(B::Zero, C::AbstractScalar) = -C

# B::AbstractScalar, x::Real
# x::Real, B::AbstractScalar
-(B::AbstractScalar, x::Real) = Scalar{typeof(value(B))}(value(B) - x)
-(x::Real, B::AbstractScalar) = Scalar{typeof(value(B))}(x - value(B))

# --- Operations involving a Scalar instance

# B::Scalar, C::Scalar
-(B::Scalar, C::Scalar) = Scalar{typeof(value(B))}(value(B) - value(C))

# --- Operations involving a One instance

# B::One, C::One
-(B::One, C::One) = Zero{typeof(value(B))}()

# B::One, C::Zero
# B::Zero, C::One
-(B::One, C::Zero) = B
-(B::Zero, C::One) = -C

# --- Operations involving a Zero instance

# B::Zero, C::Zero
-(B::Zero, C::Zero) = B
