# Implementation Details

## Blade Representation

* Multiplicative representation from Fontijne. Blades are stored as

  * a matrix containing an orthonormal basis for the subspace represented by the blade and

  * a scalar value representing the magnitude and orientation of the blade

* Using this representation, a blade's dimension ``n`` and grade ``k`` are limited by only
  the availability of RAM to store ``n \times k`` matrices.

* This representation makes it possible to leverage state-of-the-art numerical linear
  algebra algorithms to implement geometric algebra operations.
