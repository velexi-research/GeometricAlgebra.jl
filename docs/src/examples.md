# Examples

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

* Compute the outer product between blades.

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

  ```julia
  julia> A = [0, 0, 0, 2];

  julia> C = A ∧ B
  Blade{Float64}(4, 3, [-0.6 0.0 0.0; 0.8 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0], -10.0)

  julia> grade(Blade(A))
  1

  julia> grade(B)
  2

  julia> grade(C)
  3
  ```

* Compute the additive and multiplicative inverses of a vector (1-blade).

  ```julia
  julia> v = Blade([3, 4, 0, 0])
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;;], 5.0)

  julia> -v
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;;], -5.0)

  julia> inv(v)
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;;], 0.2)

  julia> 1 / v
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;;], 0.2)

  julia> v^-1
  Blade{Float64}(4, 1, [0.6; 0.8; 0.0; 0.0;;], 0.2)
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
  Blade{Float64}(3, 1, [0.8; 0.6; 0.0;;], -5.0)

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
