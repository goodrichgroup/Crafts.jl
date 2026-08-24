# Stability

Linearising the rate equations about equilibrium gives the modes the system relaxes through and the rates it relaxes at.
All of this requires the network to be detailed balanced, and throws if it is not.

## The stability matrix

[`stabilitymatrix`](@ref) is the derivative of the rate equations with respect to the densities, evaluated at equilibrium:

```jldoctest stab
julia> using Crafts, Roly, LinearAlgebra

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike(); embed3d=true));

julia> net = ReactionNetwork(asys);

julia> ξ = topotentials(asys, fill(0.1, 3), fill(8.0, 2));

julia> S = stabilitymatrix(net, ξ);

julia> size(S)
(6, 6)
```

Its eigenvalues are the relaxation rates.

Two bases are available.
[`DensityBasis`](@ref) is the default and works in the densities themselves.
[`SymmetricBasis`](@ref) rescales by the square root of the equilibrium densities, which makes the matrix symmetric and its eigendecomposition real and orthogonal:

```jldoctest stab
julia> Ssym = Symmetric(stabilitymatrix(net, ξ, SymmetricBasis()));

julia> issymmetric(Ssym)
true
```

Use that basis to look at the modes, since the density basis is not symmetric and its eigenvalues come back complex.
Conserved quantities, one per species, show up as zero eigenvalues:

```jldoctest stab
julia> count(<(-1e-8), eigvals(Ssym)) == nstructures(asys) - nspecies(asys)
true
```

The two are related by a similarity transform with `D = diagm(sqrt.(densities(asys, ξ)))`, so they share eigenvalues, and an eigenvector `v` of the symmetric form maps to `D * v`.

In-place versions [`stabilitymatrix!`](@ref) and [`stabilityjacobian!`](@ref) write into an array you supply.
A `scale` keyword multiplies the result.

## Correlation time

[`correlationtime`](@ref) is the inverse of the slowest relaxation rate, which is the timescale the system takes to equilibrate:

```jldoctest stab
julia> round(correlationtime(net, ξ); digits=4)
9.4087
```

It counts one conserved quantity per species by default; pass `nconserved` if that is not right.
It can also be applied to a matrix directly, in which case wrap it in `Symmetric` if it is one, since that route is both faster and more accurate.
An equilibrium with no decaying modes returns `Inf`.

## Derivatives

[`stabilityjacobian`](@ref) differentiates the stability matrix with respect to `ξ`:

```jldoctest stab
julia> size(stabilityjacobian(net, ξ))
(6, 6, 5)
```
