"""
reject.jl defines methods for the reject(x, y) function

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Exports

export reject

# --- Method definitions

using LinearAlgebra: ⋅

"""
    reject(vectors::Matrix, B::AbstractBlade; normalize=false)::Matrix

Compute rejections of `vectors` from `B`. When `normalize` is true, the
rejection vectors are normalized.
"""
function reject end

# ------ Specializations involving a Blade instance

# vectors::Matrix, B::Blade
function reject(vectors::Matrix, B::Blade; normalize::Bool=false)
    # --- Check arguments

    if size(vectors, 1) != dim(B)
        throw(DimensionMismatch("`dim(vectors)` not equal to `dim(B)`"))
    end

    # --- Compute rejections using modified Gram-Schmidt algorithm

    # Initialize rejections
    rejections = Matrix{typeof(volume(B))}(vectors)

    # Remove basis(B) from rejections
    for idx_B in 1:grade(B)
        B_column = basis(B)[:, idx_B]
        for idx in 1:size(rejections, 2)
            rejections[:, idx] -=
                (rejections[:, idx] ⋅ B_column) * B_column
        end
    end

    # Normalize rejection vectors
    if normalize
        for idx in 1:size(rejections, 2)
            norm = LinearAlgebra.norm(rejections[:, idx])
            if norm > 0
                rejections[:, idx] /= norm
            end
        end
    end

    return rejections
end
