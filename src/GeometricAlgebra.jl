"""
The GeometricAlgebra.jl module demonstrates a Julia module.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the XYZ package. It is subject to
the license terms in the LICENSE file found in the top-level directory of
this distribution. No part of the XYZ package, including this file, may be
copied, modified, propagated, or distributed except according to the terms
contained in the LICENSE file.
------------------------------------------------------------------------------
"""
module GeometricAlgebra

# --- Imports

import LinearAlgebra


# --- Exports

# Types
export Blade

# Functions
export norm


# --- Types

"""
    struct Blade{T<:AbstractFloat}

The Blade type represents a blade that is stored with the floating-point
precision of type `T`.
"""
struct Blade{T<:AbstractFloat}
    basis::Matrix{T}
    norm::T

    """
        Blade{T}(vectors::Array{T,2})

    Construct a Blade from a collection of vectors stored as the columns
    of a 2-dimensional array.
    """
    function Blade{T}(vectors::Array{T,2}) where {T<:AbstractFloat}
        dims = size(vectors)
        if dims[1] < dims[2]
            if dims[1] == 1
                # `vectors` is a single row vector, so convert it to a column
                # vector and call constructor for single column vector.
                return Blade{T}(reshape(vectors, dims[2]))
            else
                return new(Array{T}(undef, 0, 2), 0)
            end
        else
            F = LinearAlgebra.qr(vectors)
            basis::Matrix{T} = F.Q
            norm::T = abs(prod(LinearAlgebra.diag(F.R)))
            new(basis, norm)
        end
    end

    """
        Blade{T}(vectors::Array{T,1})

    Construct a Blade from a single vector (1-dimensional Array)
    """
    function Blade{T}(vector::Array{T,1}) where {T<:AbstractFloat}
        norm::T = LinearAlgebra.norm(vector)
        basis::Matrix{T} = reshape(vector, length(vector), 1) / norm
        new(basis, norm)
    end
end

"""
    Blade(vectors)

Construct a Blade from a collection of vectors stored as the columns of a
2-dimensional array of floating-point values.

The precision of the Blade is inferred precision from the precision of the
`vectors` Array.
"""
Blade(vectors::Array{T}) where {T<:AbstractFloat} = Blade{T}(vectors)

"""
    Blade(vectors::Array{Int})
    Blade{T}(vectors::Array{Int}) where {T<:AbstractFloat}

Construct a Blade from a collection of vectors stored as the columns of a
2-dimensional array of integer values.

When the precision of the Blade is not explicitly specified, it defaults to
Float64.
"""
Blade(vectors::Array{Int}) = Blade(convert(Array{Float64}, vectors))
Blade{T}(vectors::Array{Int}) where {T<:AbstractFloat} =
    Blade{T}(convert(Array{T}, vectors))


# --- Functions

end  # End of GeometricAlgebra.jl module
