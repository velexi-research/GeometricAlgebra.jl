"""
The operators_mgs.jl submodule defines versions of operations on subtypes of
AbstractBlade based on modified Gram-Schmidt orthogonalization (as opposed to
Householder triangularization).

Notes
-----
* Compared with algorithms based on Householder triangularization, algorithms
  based on modified Gram-Schmidt orthogonalization typically have higher
  computational performance at the cost of lower accuracy and higher memory
  usage.

------------------------------------------------------------------------------
COPYRIGHT/LICENSE. This file is part of the GeometricAlgebra.jl package. It
is subject to the license terms in the LICENSE file found in the top-level
directory of this distribution. No part of the GeometricAlgebra.jl package,
including this file, may be copied, modified, propagated, or distributed
except according to the terms contained in the LICENSE file.
------------------------------------------------------------------------------
"""
# --- Imports

# Standard library
import LinearAlgebra


# --- Binary operations

# Exports
export dual_mgs, outer_mgs

"""
    outer_mgs(B::Blade, C::Blade)

Return the outer product of `B` and `C`.

Notes
-----
* The extension of the basis of `B` to a basis for span(B, C) is computed
  using modified Gram-Schmidt orthogonalization.
"""
function outer_mgs(B::Blade, C::Blade)
    # Compute the rejections of the basis of `C` from the subspace represented
    # by `B`
    rejections = rejection(basis(C), B)

    # Orthogonalize the rejections using the modified Gram-Schmidt algorithm
    for j in 1:size(rejections, 2)
        norm = LinearAlgebra.norm(rejections[:, j])
        if norm > 0
            rejections[:, j] /= norm
        else
            return zero(B)
        end

        column = rejections[:, j]
        for j_inner in (j + 1):size(rejections, 2)
            rejections[:, j_inner] -=
                (rejections[:, j_inner] ⋅ column) * column
        end
    end

    # Construct the Blade representing the outer product
    B_wedge_C_basis = hcat(basis(B), rejections)
    Blade{typeof(B.volume)}(dim(B), grade(B) + grade(C),
                            B_wedge_C_basis, volume(B) * volume(C),
                            copy_basis=false)

#    rejections = rejection(basis(C), B)
#    Blade(hcat(basis(B), rejections), volume=volume(B) * volume(C))
end

"""
    dual_mgs(B::AbstractBlade, C::AbstractBlade)

Return the dual `B` relative to `C`.

Notes
-----
* `dual(B, C)` is only defined if (1) `B` and `C` are extended from real
  vector spaces of the same dimension and (2) the subspace represented by `B`
  is contained in subspace represented by `C`.

* The volume of `C` is ignored.

* The extension of `basis(B)` to `basis(C)` is computed using modified
  Gram-Schmidt orthogonalization.
"""
function dual_mgs(B::Blade, C::Blade)
    # --- Handle edge cases

    # Check that B and C are extended from the real vector spaces of the same
    # dimension
    if dim(B) != dim(C)
        throw(DimensionMismatch("`dim(B)` not equal to `dim(C)`"))
    end

    # Check that B is contained in C
    projection_coefficients = transpose(basis(C)) * basis(B)
    if LinearAlgebra.norm(projection_coefficients)^2 ≉ grade(B)
        throw(DimensionMismatch("`B` not contained in `C`"))
    end

    # Subspaces represented by B and C are the same
    if grade(B) == grade(C)
        return Scalar(volume(B))
    end

    # --- Compute dual using the basis vectors of `C` with the smallest
    #     projections (i.e., largest rejections) onto the basis of `B`.

    permutation = sortperm(sum(abs.(projection_coefficients), dims=2)[:, 1])
    dual_basis = rejection(basis(C)[:, permutation[1:grade(C) - grade(B)]], B)
    Blade(dual_basis, volume=volume(B))
end
