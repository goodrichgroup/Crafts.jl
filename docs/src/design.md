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

## Maximizing yield within a bond energy budget

[`maxyielddesign`](@ref) spends a budget of bond energy to make the target as abundant as possible.

```jldoctest design
julia> using Crafts, Roly, Convex, Clarabel

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false);

julia> ξ, residual = maxyielddesign(asys, 6; energy_budget=8.0, optimizer=Clarabel.Optimizer);

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
julia> ξ, ε̄ = minenergydesign(asys, 6; minyield=0.9, optimizer=Clarabel.Optimizer);

julia> round(ε̄; digits=3)
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
Here `energy_measure` picks what is minimized rather than what is capped, and the returned `ε̄` is that measure at the optimum.

## Designability

[`lineardesign`](@ref) asks a coarser question: is there any choice of parameters at all that isolates the target?
It places the target at zero excess free energy and pushes everything else below.
When that is impossible, the target set is not *designable*, and with `refine=true` the function reports the smallest designable set containing it.

## What is designable at all

Before optimizing, it is worth asking what the binding rules permit.
Pushing the parameters to their limit suppresses every structure whose composition satisfies `Mₛ·ξ < 0`, so the directions worth taking are those in the *constraint cone* `{ξ : Mₛ·ξ ≤ 0 for all s}`.
A direction on the boundary of that cone assembles exactly the structures whose constraint is tight there, and the faces of the cone therefore enumerate what can be assembled together.

Take a threefold-symmetric triangle capped by up to three copies of a second species, with a single shared bond energy:

```jldoctest cone
julia> using Crafts, Roly

julia> central = PolygonParticleSpecies(3; labels=[1, 1, 1]);

julia> outer = PolygonParticleSpecies(3);

julia> rules = BindingRules([1 1 2 1; 1 2 2 1; 1 3 2 1], [central, outer]);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false);

julia> nparticles.(structures(asys))
5-element Vector{Int64}:
 1
 1
 2
 3
 4
```

The five structures are the two monomers, and the central triangle carrying one, two or three caps.
`bondgroups=:uniform` ties the three bond types to a single energy, leaving three parameters: two chemical potentials and one `ε`.

```jldoctest cone
julia> cone = constraintcone(asys; bondgroups=:uniform)
ConstraintCone[d=3, 5 structures, 3 designable, 3 rays]

julia> designablestructures(asys; bondgroups=:uniform)
3-element Vector{Int64}:
 1
 2
 5
```

Only three of the five constraints bound the cone, so only those structures can be assembled on their own.
The dimer and trimer cannot: their constraints touch the cone without bounding it.

```jldoctest cone
julia> isdesignable(cone, 3)
false

julia> minimaldesignableset(cone, 3)
4-element Vector{Int64}:
 1
 3
 4
 5

julia> necessarychimeras(cone, 3)
3-element Vector{Int64}:
 1
 4
 5
```

Assembling the dimer means tolerating the central monomer, the trimer and the tetramer.
Designable sets are closed under intersection, so the smallest one containing a given target is unique, and [`necessarychimeras`](@ref) reports what it forces on you.

[`designablesets`](@ref) lists them all. Small sets lie on the cheap end of the face lattice, so `maxsize` keeps this affordable when the full lattice is not:

```jldoctest cone
julia> designablesets(cone; maxsize=2)
5-element Vector{Vector{Int64}}:
 [1]
 [2]
 [5]
 [1, 2]
 [2, 5]
```

## Relative yields within a designable set

Aligning the parameters with a face fixes which structures survive, but not their proportions.
Those are governed by the directions perpendicular to the face, and there are usually very few.

```jldoctest cone
julia> set = [1, 3, 4, 5];

julia> codimension(asys, set; bondgroups=:uniform)
2

julia> relativeyielddofs(asys, set; bondgroups=:uniform)
1
```

The face has codimension two, but one of those directions only rescales the whole set, leaving a single knob.
[`relativeyields`](@ref) traces what that knob can reach — every combination below is thermodynamically allowed, and nothing else is:

```jldoctest cone
julia> round.(relativeyields(asys, set, [0.0]; bondgroups=:uniform); digits=4)
4-element Vector{Float64}:
 0.125
 0.375
 0.375
 0.125

julia> round.(relativeyields(asys, set, [3.0]; bondgroups=:uniform); digits=4)
4-element Vector{Float64}:
 0.0
 0.0006
 0.0413
 0.9581
```

At the origin the proportions are the bare partition functions, here the binomial weights for occupying three equivalent sites.
Increasing the coordinate moves the population towards the tetramer.

This prediction is independent of the entropy model's absolute scale and of how far the parameters are from their asymptotic limit, which makes it a robust thing to compare against experiment.

## Polyhedral primitives

The cone functions above are built on general polyhedral computation, exposed directly for other uses.
These work in exact rational arithmetic with [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html) and need no solver.

[`extremerays`](@ref) converts inequalities into generators, and [`facetsof`](@ref) goes the other way, also reporting which generators lie on which facet.

```jldoctest polyhedra
julia> using Crafts

julia> A = Rational{BigInt}[-1 0; 0 -1; 1 1];

julia> b = Rational{BigInt}[0, 0, 1];

julia> rays, verts = extremerays(A, b);

julia> verts
3×2 Matrix{Rational{BigInt}}:
 1  0
 0  1
 0  0

julia> _, _, incidence = facetsof(verts; rays=false);

julia> Int.(incidence)
3×3 Matrix{Int64}:
 1  1  0
 0  1  1
 1  0  1
```

[`removeredundancy`](@ref) drops inequalities implied by the others, and reports which rows survived.

A cone can have far more generators than fit comfortably in memory.
[`foreachray`](@ref) hands them to a function one at a time instead of collecting them, which keeps a run that would otherwise need tens of gigabytes within a couple:

```julia
count = Ref(0)
foreachray(A, b) do generator, isray
    count[] += 1
end
```

### Faces

[`facelattice`](@ref) enumerates the faces of the cone from that incidence matrix.

```jldoctest polyhedra
julia> L = facelattice(incidence)
FaceLattice[f=1,3,3,1]

julia> fvector(L)
4-element Vector{Int64}:
 1
 3
 3
 1
```

The f-vector counts faces by dimension, starting at the empty face: a triangle has three vertices, three edges, and itself.
[`faces`](@ref) returns the faces of one dimension, each a `BitVector` over the generators.

The number of faces grows exponentially with dimension, and design cones are strongly degenerate, so the whole lattice is often out of reach.
`maxdim` stops the enumeration early and is usually what you want:

```jldoctest polyhedra
julia> fvector(facelattice(incidence; maxdim=1))
2-element Vector{Int64}:
 1
 3
```

Set `dual=true` to work down from the facets instead of up from the generators, `covers=true` to also get the cover relations, and `facetsets=true` to tag each face with the facets containing it.
