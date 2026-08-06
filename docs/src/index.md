# Crafts.jl

Crafts computes what a set of self-assembling particles does: which structures form, in what proportion at equilibrium, how fast they get there, and how stable that equilibrium is.

The particles and the rules for how they stick together come from [Roly.jl](https://github.com/goodrichgroup/Roly.jl).
Crafts takes those rules, enumerates every structure they allow, and gives each one a partition function.
From there you can ask for equilibrium yields, simulate the assembly kinetics, or look at the relaxation modes.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/goodrichgroup/Roly.jl")
Pkg.add(url="https://github.com/goodrichgroup/Crafts.jl")
```

## A first example

Three species of square particle that bind into a chain:

```jldoctest quickstart
julia> using Crafts, Roly

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()))
AssemblySystem[n=3, k=2]
```

`n` is the number of species and `k` the number of bond types.
The structures those rules allow are enumerated on first use:

```jldoctest quickstart
julia> nstructures(asys)
6

julia> nparticles.(structures(asys))
6-element Vector{Int64}:
 1
 1
 1
 2
 2
 3
```

Equilibrium yields follow from the particle concentrations and the bond energies:

```jldoctest quickstart
julia> ϕs, εs = fill(0.1, 3), fill(8.0, 2);

julia> round.(yields(asys, ϕs, εs); digits=4)
6-element Vector{Float64}:
 0.1063
 0.0144
 0.1063
 0.092
 0.092
 0.5891
```

The last entry is the full three-particle chain, so 59% of the structures present are the target.

To follow the assembly in time instead, build a reaction network and integrate it:

```jldoctest quickstart
julia> net = ReactionNetwork(asys)
ReactionNetwork[4/4 reactions active]

julia> ts, us = simulate_kinetics(net, ϕs, εs; T=1000.0);

julia> round.(us[:, end]; digits=4)
6-element Vector{Float64}:
 0.0135
 0.0018
 0.0135
 0.0117
 0.0117
 0.0748
```

Those are number densities rather than yields, and they are the equilibrium the yields above describe.

## Where to go next

  - [Assembly systems](assemblysystems.md) — enumerating structures and their partition functions.
  - [Bond potentials](potentials.md) — modelling the interaction, and the solvers that integrate it out.
  - [Equilibrium](equilibrium.md) — yields, densities, and their derivatives.
  - [Kinetics](kinetics.md) — reaction networks, rate kernels, and time integration.
  - [Stability](stability.md) — relaxation modes around equilibrium.
  - [Diffusion and rates](rates.md) — hydrodynamics and diffusion-limited binding rates.
