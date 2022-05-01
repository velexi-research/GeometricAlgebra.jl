# GeometricAlgebra.jl

[GeometricAlgebra.jl](https://github.com/velexi-corporation/GeometricAlgebra.jl)
defines a collection of basic types and operations that support numerical
geometric algebra computations. Our aim is to simplify the process of implementing
numerical algorithms for geometric operations expressed algebraically in the language of
geometric algebra.

The GeometricAlgebra package features:

* an easy-to-use interface for constructing geometric algebra objects,

* intuitive notation for operations on and between geometric algebra objects, and

* behavior that is consistent with expectations for numerical linear algebra libraries.

## Overview

Geometric algebra extends classical linear algebra (of inner product spaces over the real
numbers ``\mathbb{R}``) by including the higher-dimensional siblings of the vector (a
1-dimensional object) and defining algebraic operations that augment the basic operations
on inner product spaces, such as the exterior product and geometric product. With these
extensions, geometric algebra makes it possible to use simple algebraic expressions to
represent geometric objects (e.g., subspaces) and operations (e.g., projection, computing
the orthogonal complement).

The key mathematical object introduced in geometric algebra is the _blade_. Geometrically,
blades represent two concepts (1) ``k``-dimensional subspaces of ``\mathbb{R}^n`` and
(2) ``k``-dimensional hypervolumes (of arbitrary shape). The _exterior product_
(``\wedge``), which combines two blades and yields a blade, is an example of an algebraic
operation with geometric meaning. For instance, it enables any blade ``B`` to be expressed
_algebraically_ in terms of a basis ``b_1, \ldots b_k`` for the subspace represented by
the blade:

```math
B = b_1 \wedge \cdots \wedge b_k.
```

The GeometricAlgebra package implements types and functions that make it straightforward
to do _computational_ geometric algebra calculations.

## Getting Started

* Add the Velexi Julia package registry.

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> registry add https://github.com/velexi-corporation/JuliaRegistry.git
  ```

  !!! tip "Only needed once"

      This step only needs to be performed once per Julia installation.

  !!! note "GeometricAlgebra is registered with a local Julia package registry"

      The Velexi registry needs to be added to your Julia installation because
      GeometricAlgebra is currently registered with Velexi Julia package registry (not
      the General Julia package registry).

* Install the GeometricAlgebra package via the Pkg REPL. That's it!

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> add GeometricAlgebra
  ```

## Summary of Types

GeometricAlgebra defines concrete types for the following hierarchy of mathematical objects.

* multivector: [`Multivector`](@ref)
  * blade: [`Blade`](@ref)
    * scalar (``0``-blade) : [`Scalar`](@ref)
      * multiplicative identity: [`One`](@ref)
      * additive identity: [`Zero`](@ref)
    * pseudoscalar (``n``-blade): [`Pseudoscalar`](@ref)

## Summary of Operations

### Unary Operations

* additive inverse: [`-`](@ref -(::AbstractMultivector))
* multiplicative inverse: [`inv`](@ref)
* reverse: [`reverse`](@ref)
* dual: [`dual`](@ref)

### Binary Operations

* addition/subtraction: [`+`](@ref Base.:(+)), [`-`](@ref)
* multiplication/division (geometric product): [`*`](@ref), [`/`](@ref)
* exterior product: [`∧`](@ref) ([`wedge`](@ref))
* inner product (left contraction):
  [`<`](@ref GeometricAlgebra.contract_left) ([`contract_left`](@ref)),
  [`⋅`](@ref GeometricAlgebra.contract_left) ([`dot`](@ref GeometricAlgebra.contract_left))

### Geometric Operations

* projection: [`project`](@ref)

## Supported Numeric Types

All GeometricAlgebra types are parametrized by the numeric type used to represent the real
numbers. Only subtypes of `AbstractFloat` are supported (i.e., `Float16`, `Float32`,
`Float64`, `BigFloat`). Other subtypes of `Real` (e.g., `Rational`, `Irrational`) are not
allowed as parameters to GeometricAlgebra types.

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
