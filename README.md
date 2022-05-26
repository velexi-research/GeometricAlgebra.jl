# GeometricAlgebra.jl

[------------------------------------ BADGES: BEGIN ------------------------------------]: #

<table>
  <tr>
    <td>Documentation</td>
    <td>
      <a href="https://velexi-research.github.io/GeometricAlgebra.jl/dev/"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/docs-dev-blue.svg"/></a>
      <a href="https://velexi-research.github.io/GeometricAlgebra.jl/stable/"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/docs-stable-blue.svg"/></a>
    </td>
  </tr>

  <tr>
    <td>Build Status</td>
    <td>
      <a href="https://github.com/velexi-research/GeometricAlgebra.jl/actions/workflows/CI.yml"><img style="vertical-align: bottom;" src="https://github.com/velexi-research/GeometricAlgebra.jl/actions/workflows/CI.yml/badge.svg"/></a>
      <a href="https://codecov.io/gh/velexi-research/GeometricAlgebra.jl">
        <img style="vertical-align: bottom;" src="https://codecov.io/gh/velexi-research/GeometricAlgebra.jl/branch/main/graph/badge.svg?token=2PSLG9EGAK"/>
      </a>
    </td>
  </tr>

  <!-- Miscellaneous Badges -->
  <tr>
    <td colspan=2 align="center">
      <a href="https://github.com/velexi-research/GeometricAlgebra.jl/issues"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat"/></a>
      <a href="https://github.com/invenia/BlueStyle"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/code%20style-blue-4495d1.svg"/></a>
      <a href="http://hits.dwyl.com/velexi-research/GeometricAlgebrajl"><img style="vertical-align: bottom;" src="https://hits.dwyl.com/velexi-research/GeometricAlgebrajl.svg?style=flat-square&show=unique"/></a>
    </td>
  </tr>
</table>

[------------------------------------- BADGES: END -------------------------------------]: #

GeometricAlgebra defines a collection of basic types and operations that support numerical
geometric algebra computations. Our aim is to simplify the process of implementing
numerical algorithms for geometric operations expressed algebraically in the language of
geometric algebra.

The GeometricAlgebra package features:

* an easy-to-use interface for constructing geometric algebra objects,

* intuitive notation for operations on and between geometric algebra objects, and

* behavior (including limitations) analogous to standard numerical linear algebra
  libraries.

## Overview

Geometric algebra extends classical linear algebra (of inner product spaces over the real
numbers $\mathbb{R}$) by including the higher-dimensional siblings of the vector (a
1-dimensional object) and defining algebraic operations that augment the basic operations
on inner product spaces (i.e., vector addition, scalar multiplication, and inner product).
With these extensions, geometric algebra makes it possible to use simple algebraic
expressions to represent geometric objects (e.g., subspaces) and operations (e.g.,
projection, computing the orthogonal complement).

The key mathematical object introduced in geometric algebra is the _blade_. Geometrically,
blades represent two concepts (1) subspaces of $\mathbb{R}^n$ and (2) $k$-dimensional
hypervolumes (with dimension less than or equal to $n$ and arbitrary shape). The
_exterior product_ ($\wedge$), which combines lower-grade blades into higher-grade blades
(or zero), is an example of an algebraic operation with geometric meaning. For instance,
it enables any blade $B$ to be expressed _algebraically_ in terms of a basis
$b_1, \ldots b_k$ for the subspace represented by the blade:

\[
B = b_1 \wedge \cdots \wedge b_k.
\]

The GeometricAlgebra package implements types and functions that make it straightforward
to do _computational_ geometric algebra calculations.

## Getting Started

* Add the Velexi Julia package registry.

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> registry add https://github.com/velexi-research/JuliaRegistry.git
  ```

  __Notes__

  * _Only needed once_. This step only needs to be performed once per Julia installation.

  * _GeometricAlgebra is registered with a local Julia package registry_.
    The Velexi registry needs to be added to your Julia installation because
    GeometricAlgebra is currently registered with Velexi Julia package registry (not
    the General Julia package registry).

* Install the GeometricAlgebra package via the Pkg REPL. That's it!

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> add GeometricAlgebra
  ```

## Examples

* Create a blade.

  ```julia
  julia> vectors = Matrix([3 3; -4 -4; 0 1]);

  julia> B = Blade(vectors)
  Blade{Float64}(3, 2, [-0.6 0.0; 0.8 0.0; 0.0 -1.0], 5.0)

  julia> dim(B)
  3

  julia> grade(B)
  2

  julia> norm(B)
  5.0
  ```

* Compute the outer product of vectors (1-blades).

  ```julia
  julia> u = [3, -4, 0, 0];  # a vector is a 1-blade

  julia> v = [3, -4, 1, 0];  # a vector is a 1-blade

  julia> B = u ∧ v
  Blade{Float64}(3, 2, [-0.6 0.0; 0.8 0.0; 0.0 -1.0], 5.0)

  julia> grade(Blade(u))
  1

  julia> grade(Blade(v))
  1

  julia> grade(B)
  2
  ```

* Compute the outer product of blades.

  ```julia
  julia> A = Blade([0, 0, 0, 2]);

  julia> B = Blade(Matrix([3 3; -4 -4; 0 1; 0 0]));

  julia> C = A ∧ B
  Blade{Float64}(4, 3, [0.0 0.6 0.0; 0.0 -0.8 0.0; 0.0 0.0 1.0; -1.0 0.0 0.0], -10.0)

  julia> grade(A)
  1

  julia> grade(B)
  2

  julia> grade(C)
  3
  ```

* Compute the additive and multiplicative inverses of a vector (1-blade).

  ```julia
  julia> v = Blade([3, 4, 0, 0])
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;], 5.0)

  julia> -v
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;], -5.0)

  julia> inv(v)
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;], 0.2)

  julia> 1 / v
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;], 0.2)

  julia> v^-1
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;], 0.2)
  ```

* Create a pseudoscalar for $\mathbb{G}^n$.

  ```julia
  julia> n = 10;
  julia> value = -2;
  julia> B = Pseudoscalar(n, value)
  Pseudoscalar{Float64}(10, -2.0)

  julia> grade(B)
  10

  julia> norm(B)
  2.0

  julia> volume(B)
  -2.0
  ```

* Compute the dual of a blade.

  ```julia
  julia> vectors = Matrix([3 3; -4 -4; 0 1]);

  julia> B = Blade(vectors);

  julia> dual(B)
  Blade{Float64}(3, 1, [0.8; 0.6; 0.0;], -5.0)

  julia> isapprox(dot(B.basis[:, 1], dual(B).basis), 0; atol=1e-15)
  true

  julia> isapprox(dot(B.basis[:, 2], dual(B).basis), 0; atol=1e-15)
  true

  julia> B ∧ dual(B)
  Pseudoscalar{Float64}(3, 25.0)
  ```

  We can confirm that is the dual of blade $B$ is the same as division by the unit
  pseudoscalar for $\mathbb{G}^n$.

  ```julia
  julia> dual(B) ≈ B / Pseudoscalar(3, 1)
  true
  ```

* Compute the geometric product of two vectors. Note that the geometric product of two
  vectors is a _multivector_ (a linear combination of blades).

  ```julia
  julia> u = Blade([1, 0, 0]);

  julia> v = Blade([sqrt(3), 1, 0]);

  julia> M = u * v
  Multivector{Float64}(3, DataStructures.SortedDict{Int64, Vector{AbstractBlade}, Base.Order.ForwardOrdering}(0 => [Scalar{Float64}(1.7320508075688776)], 2 => [Blade{Float64}(3, 2, [1.0 0.0; 0.0 1.0; 0.0 0.0], 1.0)]), 2.8284271247461903)

  julia> M[0]
  1-element Vector{AbstractBlade}:
   Scalar{Float64}(1.7320508075688776)

  julia> M[1]
  AbstractBlade[]

  julia> M[3]
  1-element Vector{AbstractBlade}:
   Blade{Float64}(3, 2, [1.0 0.0; 0.0 1.0; 0.0 0.0], 1.0)
  ```

  Right multiplication of `M` by the multiplicative inverse of `v` yields `u`.

  ```julia
  julia> M * inv(v) ≈ u
  true

  julia> M * (1 / v) ≈ u
  true
  ```

  Left multiplication of `M` by the multiplicative inverse of `u` yields `v`.

  ```julia
  julia> inv(u) * M ≈ v
  true

  julia> (1 / u) * M ≈ v
  true
  ```

## Related Packages

There are several active geometric algebra packages in the Julia ecosystem. Most of the
available packages are implemented using an _additive_ blade representation and emphasize
_symbolic_ computations. To the best of our knowledge, GeometricAlgebra.jl is the only
package currently uses a _multiplicative_ blade representation and focuses on _numerical_
computations.

* [Grassmann.jl](https://grassmann.crucialflow.com/)

* [Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl)

* [Jollywatt/GeometricAlgebra.jl](https://github.com/Jollywatt/GeometricAlgebra.jl)

* [serenity4/GeometricAlgebra.jl](https://github.com/serenity4/GeometricAlgebra.jl)

* [GAlgebra.jl](https://github.com/pygae/GAlgebra.jl)
