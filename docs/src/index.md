# Crafts.jl

Crafts.jl (**C**ombinato**r**ial **A**nalysis **F**ramework for **T**argeted **S**elf-assembly) computes and optimizes the equilibrium properties and assembly kinetics of mixtures of self-assembling particles.

Starting from a set of binding rules defined within its sister package [Roly.jl](https://github.com/goodrichgroup/Roly.jl), Crafts.jl enumerates every structure those rules allow and then predicts their number densities and yields.
It further allows investigation of the assembly kinetics by integrating the rate equations governing the assembly.

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

Equilibrium number densities follow from the particle concentrations and the bond energies:

```jldoctest quickstart
julia> ϕs, εs = fill(0.1, 3), fill(8.0, 2);

julia> round.(densities(asys, ϕs, εs); digits=4)
6-element Vector{Float64}:
 0.0135
 0.0018
 0.0135
 0.0117
 0.0117
 0.0748
```

The last entry is the full three-particle chain.
[`yields`](@ref) gives the same thing as fractions of all structures present, which puts the chain at 59%.

To follow the assembly in time instead, build a reaction network and integrate the rate equations:

```jldoctest quickstart
julia> net = ReactionNetwork(asys)
ReactionNetwork[4/4 reactions active]

julia> ts, ρs = simulate_kinetics(net, ϕs, εs; T=1000.0);

julia> round.(ρs[:, end]; digits=4)
6-element Vector{Float64}:
 0.0135
 0.0018
 0.0135
 0.0117
 0.0117
 0.0748
```

`ρs` holds number densities, one column per time point in `ts`, so the last column is the equilibrium computed above.

## Where to go next

  - [Assembly systems](assemblysystems.md) — enumerating structures and their partition functions.
  - [Bond potentials](potentials.md) — modelling the interaction, and the solvers that integrate it out.
  - [Equilibrium](equilibrium.md) — yields, densities, and their derivatives.
  - [Kinetics](kinetics.md) — reaction networks, rate kernels, and time integration.
  - [Stability](stability.md) — relaxation modes around equilibrium.
  - [Diffusion and rates](rates.md) — hydrodynamics and diffusion-limited binding rates.
