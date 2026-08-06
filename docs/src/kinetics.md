# Kinetics

Equilibrium says where the system ends up.
To see how it gets there, build a reaction network and integrate it.

## Reaction networks

A [`ReactionNetwork`](@ref) enumerates every way one structure can split into two, and so every way two can combine into one:

```jldoctest kin
julia> using Crafts, Roly

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()));

julia> net = ReactionNetwork(asys)
ReactionNetwork[4/4 reactions active]

julia> reactions(net)
4-element Vector{Tuple{Int64, Int64, Int64}}:
 (1, 2, 4)
 (2, 3, 5)
 (1, 5, 6)
 (4, 3, 6)
```

Each triple `(i, j, k)` means structures `i` and `j` combine into `k`, indexing into `structures(asys)`.
Pass `maxbonds` to ignore reactions that break more bonds than that at once.

Iterating the network gives a [`Reaction`](@ref) per entry, which knows more about itself:

```jldoctest kin
julia> rxn = net[3]
Reaction[1 + 5 ⟷ 6, 1 bond(s)]

julia> indices(rxn)
(1, 5, 6)

julia> bondcounts(rxn)
2-element Vector{Int64}:
 1
 0

julia> nparticles(product(rxn))
3
```

[`reactants`](@ref) gives the two fragments as free structures, [`product`](@ref) the structure that forms, [`cut`](@ref) the bonds broken and [`halves`](@ref) which part of the product each fragment is.

## Rates

Every reaction has a forward and a backward rate, both `1.0` by default.
[`rate!`](@ref) sets them from a kernel, which is any function taking a `Reaction` and returning a rate:

```jldoctest kin
julia> rate!(net; fwdkernel=rxn -> 1.0 + nparticles(product(rxn)));

julia> fwdrates(net)
4-element Vector{Float64}:
 3.0
 3.0
 4.0
 4.0
```

The backward kernel defaults to the forward one.
Rates are attempt frequencies, so equal forward and backward kernels leave the equilibrium unchanged and only alter how fast it is reached; [`isdetailedbalanced`](@ref) reports whether that holds:

```jldoctest kin
julia> isdetailedbalanced(net)
true
```

A rate of zero deactivates a reaction.
`net.active` marks which reactions are live, and rates must be finite and non-negative.

The same kernels can be passed to the constructor: `ReactionNetwork(asys; fwdkernel, bwdkernel)`.
See [Diffusion and rates](rates.md) for physical kernels.

## Time integration

[`simulate_kinetics`](@ref) integrates the rate equations, starting from pure monomers unless told otherwise:

```jldoctest kin
julia> ts, ρs = simulate_kinetics(net, fill(0.1, 3), fill(8.0, 2); T=1000.0);

julia> size(ρs, 1) == nstructures(asys)
true

julia> round.(ρs[:, end]; digits=4)
6-element Vector{Float64}:
 0.0135
 0.0018
 0.0135
 0.0117
 0.0117
 0.0748
```

`ρs` holds number densities, one column per time point in `ts`.
Pass `saveat` for specific times, `initial_densities` for a different starting state, and any other keyword through to the ODE solver:

```jldoctest kin
julia> ts, ρs = simulate_kinetics(net, fill(0.1, 3), fill(8.0, 2); T=1000.0, saveat=0:100:1000);

julia> length(ts)
11
```
