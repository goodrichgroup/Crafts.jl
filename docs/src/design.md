# Design

Given an assembly system, the design functions find the particle concentrations and bond energies that make a chosen structure dominate the yield.
They are the inverse of the equilibrium calculation: instead of taking parameters and returning yields, they take a target and return parameters.

Both problems are convex, and they are solved through [Convex.jl](https://jump.dev/Convex.jl/stable/).
Crafts does not bundle a solver, so load Convex together with one of your choice.
[Clarabel](https://clarabel.org) is a good default: it is written in Julia, so it adds no binary dependency.

```julia
using Pkg
Pkg.add(["Convex", "Clarabel"])
```

The design functions are unavailable until `Convex` is loaded, and every one of them takes an `optimizer`.
[`maxyielddesign`](@ref) and [`minenergydesign`](@ref) additionally require `maxdensity`, the total particle density to stay below.

## Maximizing yield within a bond energy budget

[`maxyielddesign`](@ref) spends a budget of bond energy to make the target as abundant as possible.

```jldoctest design
julia> using Crafts, Roly, Convex, Clarabel

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false);

julia> ξ, residual = maxyielddesign(asys, 6; maxdensity=1, energy_budget=8.0, optimizer=Clarabel.Optimizer);

julia> round.(ξ; digits=3)
5-element Vector{Float64}:
 -5.485
 -8.125
 -5.485
  8.0
  8.0
```

`ξ` holds the three chemical potentials followed by the two bond energies, in the same order [`densities`](@ref) expects.
Structure 6 is the full three-particle chain.

```jldoctest design
julia> M, Ωs = compositionmatrix(asys), partitionfunctions(asys);

julia> round.(yields(ξ, M, Ωs); digits=4)
6-element Vector{Float64}:
 0.0678
 0.0048
 0.0678
 0.0598
 0.0598
 0.74
```

The chain reaches 74%, against the 59% the same bond energy gives when the concentrations are left uniform.
`residual` is `log Σ_{j≠i} ρ_j/ρ_i`, so smaller is better.

`energy_budget` bounds `energy_measure` of the bond energies, which is `mean` by default and can be any function of them — `maximum` caps the strongest bond instead of the average.
It is a cap rather than a target: bond energy also favours the larger off-target structures, so the solution is free to come in under it, as long as nothing better lies further out.
Setting `uniform_energy=true` ties every bond type to a single energy, which is what an experiment offering only one interaction strength can realize.

Passing several indices designs for all of them at once, with `relative_yields` setting the densities they are held at relative to each other.

## Minimizing bond strength at fixed yield

[`minenergydesign`](@ref) answers the opposite question — the weakest bonds that still reach a given yield.
Strong bonds are hard to realize and slow to equilibrate, so this is usually the more practical direction.

```jldoctest design
julia> ξ, residual_energy = minenergydesign(asys, 6; maxdensity=1, minyield=0.9, optimizer=Clarabel.Optimizer);

julia> round(residual_energy; digits=3)
10.171

julia> round.(yields(ξ, M, Ωs); digits=4)
6-element Vector{Float64}:
 0.0254
 0.0007
 0.0254
 0.0243
 0.0243
 0.9
```

The yield constraint is met exactly, since making the bonds any weaker would break it.
Here `energy_measure` picks what is minimized rather than what is capped, and the returned `residual_energy` is that measure at the optimum.

## Designability

[`lineardesign`](@ref) asks a coarser question: is there any choice of parameters at all that isolates the target?
It places the target at zero excess free energy and pushes everything else below.
When that is impossible, the target set is not *designable*, and with `refine=true` the function reports the smallest designable set containing it.
