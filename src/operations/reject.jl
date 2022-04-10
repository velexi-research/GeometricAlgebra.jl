#   Copyright (c) 2020-2022 Velexi Corporation
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
reject.jl defines methods for the reject(x, y) function
"""

# --- Exports

export reject

# --- Method definitions

using LinearAlgebra: ⋅

"""
    reject(vectors::Matrix, B::AbstractBlade; normalize=false)::Matrix

Compute rejections of `vectors` from `B`. When `normalize` is true, the rejection vectors
are normalized.
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
            rejections[:, idx] -= (rejections[:, idx] ⋅ B_column) * B_column
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
