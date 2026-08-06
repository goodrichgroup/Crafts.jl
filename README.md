# Crafts.jl

[![Build Status](https://github.com/goodrichgroup/Crafts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/goodrichgroup/Crafts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/goodrichgroup/Crafts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/goodrichgroup/Crafts.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://goodrichgroup.github.io/Crafts.jl/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://goodrichgroup.github.io/Crafts.jl/dev/)

Crafts.jl computes the equilibrium composition and the assembly kinetics of mixtures of self-assembling particles.

Starting from a set of binding rules defined with [Roly.jl](https://github.com/goodrichgroup/Roly.jl), it enumerates every structure those rules allow and assigns each one an internal partition function.
From these it obtains equilibrium number densities and yields as functions of the particle concentrations and the binding energies.
The same structures, together with the reversible binary reactions connecting them, define a system of rate equations, which can be integrated in time or linearised about equilibrium to give relaxation modes and timescales.

The underlying model is a dilute solution of rigid particles that bind at discrete sites, in two or three dimensions.
Different structures interact only through the bonds they form, and the partition function of each structure is that of a rigid cluster held together by a bond potential.

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
yields(asys, ϕs, εs)                  # equilibrium yield of every structure
```

To follow the assembly in time instead, build a reaction network and integrate it:

```julia
net = ReactionNetwork(asys)
ts, us = simulate_kinetics(net, ϕs, εs; T=1000.0)
```

See the [documentation](https://goodrichgroup.github.io/Crafts.jl/dev/) for bond potentials, rate kernels, stability analysis, and the full API.
