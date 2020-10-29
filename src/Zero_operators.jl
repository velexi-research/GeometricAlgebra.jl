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

dual(B::Zero; dim::Union{Integer, Nothing}=nothing) =
    error("The dual of Zero is not well-defined")

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
proj(B::Real, C::Zero) = C

# ------ dual(B, C)

dual(B::Zero, C::Zero) = dual(B)

# Special cases
dual(B::Zero, C::Blade) = dual(B)
dual(B::Zero, C::Pseudoscalar) = dual(B)
dual(B::Zero, C::Scalar) = dual(B)
dual(B::Zero, C::One) = dual(B)
dual(B::Zero, C::Real) = dual(B)
