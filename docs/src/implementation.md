# Implementation Details

## Blade Representation

In GeometricAlgebra, blades are represented using a _multiplicative representation_ where
each blade is stored as

* a matrix containing an orthonormal basis for the subspace represented by the blade and

* a scalar value representing the magnitude and orientation of the blade.

For numerical computing, this representation has the following advantages.

* A blade's dimension ``n`` and grade ``k`` are limited by only the availability of RAM to
  store ``n \times k`` matrices.

* State-of-the-art numerical linear algebra algorithms can be leveraged to implement
  geometric algebra operations.

For more details, see Fontijne's PhD Thesis in the [References](@ref References).
