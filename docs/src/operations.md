```@meta
CurrentModule = GeometricAlgebra
```

# Unary Operations
```@docs
-(M::AbstractMultivector)
inv
reverse
dual
dual(B::AbstractBlade, C::AbstractBlade)
```

# Binary Operations

```@docs
+(M::AbstractMultivector, N::AbstractMultivector)
-(M::AbstractMultivector, N::AbstractMultivector)
*(M::AbstractMultivector, N::AbstractMultivector)
/(M::AbstractMultivector, N::AbstractMultivector)
wedge
∧
contract_left
<
dot
```

## Geometric Operations

```@docs
project
```
