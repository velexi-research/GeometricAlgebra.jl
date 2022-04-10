# GeometricAlgebra.jl

[GeometricAlgebra.jl](https://github.com/velexi-corporation/GeometricAlgebra.jl)
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
numbers ``\mathbb{R}``) by including the higher-dimensional siblings of the vector (a
1-dimensional object) and defining algebraic operations that augment the basic operations
on inner product spaces (i.e., vector addition, scalar multiplication, and inner product).
With these extensions, geometric algebra makes it possible to use simple algebraic
expressions to represent geometric objects (e.g., subspaces) and operations (e.g.,
projection, computing the orthogonal complement).

The key mathematical object introduced in geometric algebra is the _blade_. Geometrically,
blades represent two concepts (1) subspaces of ``\mathbb{R}^n`` and (2) ``k``-dimensional
hypervolumes (with dimension less than or equal to ``n`` and arbitrary shape). The
_exterior product_ (``\wedge``), which combines lower-grade blades into higher-grade blades
(or zero), is an example of an algebraic operation with geometric meaning. For instance,
it enables any blade ``B`` to be expressed _algebraically_ in terms of a basis
``b_1, \ldots b_k`` for the subspace represented by the blade:

```math
B = b_1 \wedge \cdots \wedge b_k.
```

The GeometricAlgebra package implements types and functions that make it straightforward
to do _computational_ geometric algebra calculations.

## Getting Started

Install the GeometricAlgebra package via the Pkg REPL. That's it!

```julia
pkg> add GeometricAlgebra  # Press ']' to enter the Pkg REPL mode.
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

* additive inverse: [`inverse`](@ref), [`-`](@ref GeometricAlgebra.inverse)
* multiplicative inverse: [`reciprocal`](@ref)
* reverse: [`reverse`](@ref)
* dual: [`dual`](@ref)

### Binary Operations

* addition/subtraction: [`+`](@ref Base.:(+)), [`-`](@ref)
* multiplication/division (geometric product): [`*`](@ref), [`/`](@ref)
* exterior product: [`wedge`](@ref), [`\wedge`](@ref GeometricAlgebra.wedge)
* inner product (left contraction): [`contract_left`](@ref),
  [`<`](@ref GeometricAlgebra.contract_left), [`\cdot`](@ref GeometricAlgebra.contract_left)

### Geometric Operations

* projection: [`project`](@ref)

## Supported Numeric Types

All GeometricAlgebra types are parametrized by the numeric type used to represent the real
numbers. Only subtypes of `AbstractFloat` are supported (i.e., `Float16`, `Float32`,
`Float64`, `BigFloat`). Other subttypes of `Real` (e.g., `Rational`, `Irrational`) are not
allowed as parameters to GeometricAlgebra types.