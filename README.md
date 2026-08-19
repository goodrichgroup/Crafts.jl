# Crafts.jl

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://goodrichgroup.github.io/Crafts.jl/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://goodrichgroup.github.io/Crafts.jl/dev/)
[![Build Status](https://github.com/goodrichgroup/Crafts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/goodrichgroup/Crafts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/goodrichgroup/Crafts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/goodrichgroup/Crafts.jl)

Crafts.jl (_<ins>C</ins>ombinato<ins>r</ins>ial <ins>A</ins>nalysis <ins>F</ins>ramework for <ins>T</ins>argeted <ins>S</ins>elf-assembly_) computes and optimizes the equilibrium properties and assembly kinetics of mixtures of self-assembling particles.

Starting from a set of binding rules defined within its partner package [Roly.jl](https://github.com/goodrichgroup/Roly.jl), Crafts.jl enumerates every structure those rules allow and then predicts their number densities and yields.
It further allows investigation of the assembly kinetics by integrating the rate equations governing the assembly, and inverse design of the concentrations and bond energies that make a chosen structure dominate.

## Installation

To install Crafts.jl directly from your Julia REPL, first press `]` to enter Pkg mode, and then run
```
pkg> add https://github.com/goodrichgroup/Roly.jl
pkg> add https://github.com/goodrichgroup/Crafts.jl
```

## Basic usage

Three species of square particle that bind into a chain:

```julia
using Crafts, Roly

rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
asys = AssemblySystem(rules, EntropyModel(TreeLike()))

ϕs, εs = fill(0.1, 3), fill(8.0, 2)   # particle concentrations and binding energies
densities(asys, ϕs, εs)               # equilibrium number density of every structure
yields(asys, ϕs, εs)                  # fractions of all structures present
```

To follow the assembly in time instead, build a reaction network and integrate the rate equations.
The last column of `ρs` is the equilibrium computed above:

```julia
net = ReactionNetwork(asys)
ts, ρs = simulate_kinetics(net, ϕs, εs; T=1000.0)
```

To design instead of predict, load [Convex.jl](https://jump.dev/Convex.jl/stable/) with a solver of your choice.
This finds the weakest bonds that still put 90% of the population into the three-particle chain:

```julia
using Convex, Clarabel

ξ, residual_energy = minenergydesign(asys, 6; maxdensity=1, minyield=0.9, optimizer=Clarabel.Optimizer)
```

See the [documentation](https://goodrichgroup.github.io/Crafts.jl/dev/) for bond potentials, rate kernels, stability analysis, polyhedral computation, and the full API.
