"""
Zero_oeprators.jl defines operators for the Zero type

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Unary operators from the AbstractMultivector and AbstractBlade interfaces

-(B::Zero) = B

dual(B::Zero; dim::Union{Integer, Nothing}=nothing) = dual_of_zero()

# --- Binary operators from the AbstractMultivector and AbstractBlade interfaces

# ------ +(B, C)

+(B::Zero, C::Zero) = B

# ------ -(B, C)

-(B::Zero, C::Zero) = B

# ------ *(B, C)

*(B::Zero, C::Zero) = B

# ------ /(B, C)

/(B::Zero, C::Zero) = Scalar{typeof(value(B))}(NaN)

# ------ wedge(B, C)

wedge(B::Zero, C::Zero) = B

# ------ contractl(B, C)

contractl(B::Zero, C::Zero) = B

# ------ proj(B, C)

# Operations involving Reals
proj(B::Zero, C::Real) = B
proj(B::Real, C::Zero) = C

# Operations involving Vectors
proj(B::Zero, C::Vector{<:Real}) = B
proj(B::Vector{<:Real}, C::Zero) = C

# ------ dual(B, C)

# dual(B::Zero, C)
dual_of_zero() = error("The dual of Zero is not well-defined")

dual(B::Zero, C::Blade) = dual_of_zero()
dual(B::Zero, C::Pseudoscalar) = dual_of_zero()
dual(B::Zero, C::Scalar) = dual_of_zero()
dual(B::Zero, C::One) = dual_of_zero()
dual(B::Zero, C::Real) = dual_of_zero()

# dual(B, C::Zero)
dual_relative_to_zero() =
    error("The dual of anything relative to Zero is not well-defined")

dual(M::Multivector, C::Zero) = dual_relative_to_zero()
dual(B::Blade, C::Zero) = dual_relative_to_zero()
dual(B::Pseudoscalar, C::Zero) = dual_relative_to_zero()
dual(B::Scalar, C::Zero) = dual_relative_to_zero()
dual(B::One, C::Zero) = dual_relative_to_zero()
dual(B::Zero, C::Zero) = dual_relative_to_zero()
dual(B::Real, C::Zero) = dual_relative_to_zero()
