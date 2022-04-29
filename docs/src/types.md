```@meta
CurrentModule = GeometricAlgebra
```

# Types

## Type Hierarchy

The relationships between GeometricAlgebra types mirror the hierarchical relationships
between the mathematical objects they represent.

* `AbstractMultivector` represents an arbitrary multivector. All concrete GeometricAlgebra
  types are subtypes of `AbstractMultivector`.

* `AbstractBlade` represents an arbitrary blade. Any mathematical objects that is logically
  a blade (i.e., a blade, a pseudoscalar, or a scalar) is represented by a concrete
  GeometricAlgebra type that is a subtype of `AbstractBlade`.

* `AbstractScalar` represents an arbitrary scalar value. Any mathematical objects that is
  logically scalar (i.e., a scalar or a special scalar value) is represented by a concrete
  GeometricAlgebra type that is a subtype of `AbstractScalar`.

```
AbstractMultivector
├─ Multivector
└─ AbstractBlade
   ├─ Blade
   ├─ Pseudoscalar
   └─ AbstractScalar
      ├─ Scalar
      ├─ One
      └─ Zero
```

-------------------------------------------------------------------------------------------
## `AbstractMultivector`

```@docs
AbstractMultivector
```

### Interface

The `AbstractMultivector` interface is defined by the following functions.

#### Properties

```@docs
dim
grades
blades
norm
getindex
```

#### Operations

```@docs
inverse
reverse
dual
```

#### Utility Functions

```@docs
zero
iszero
one
isone
convert
```

-------------------------------------------------------------------------------------------
## `AbstractBlade` Interface

```@docs
AbstractBlade
```

### Interface

The `AbstractBlade` interface is defined by the following functions.

#### Properties 

```@docs
grade
basis
volume
sign
```

#### Operations

```@docs
dual(B::AbstractBlade, C::AbstractBlade)
inv
```

-------------------------------------------------------------------------------------------
## `AbstractScalar` Interface

```@docs
AbstractScalar
```

The `AbstractScalar` interface is defined by the following functions.

### Properties

```@docs
value
```

-------------------------------------------------------------------------------------------
## Concrete Types

### Multivector

!!! warning
    Support for the `Multivector` type is not fully implemented yet.

```@docs
Multivector
Multivector(blades::Vector{<:AbstractBlade})
Multivector(multivectors::Vector{<:AbstractMultivector})
```

### Blade

```@docs
Blade
Blade(vectors::Matrix{<:Real})
Blade(B::Blade)
Blade(x::Real)
```

### Pseudoscalar

```@docs
Pseudoscalar
Pseudoscalar(dim::Integer, value::AbstractFloat)
Pseudoscalar(B::Pseudoscalar; value::Real=value(B))
```

### Scalar

```@docs
Scalar
Scalar(value::AbstractFloat)

```

### One

```@docs
One
One()
```

### Zero

```@docs
Zero
Zero()
```

-------------------------------------------------------------------------------------------
