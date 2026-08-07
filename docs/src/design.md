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

## Maximizing yield at fixed bond strength

[`convexdesign`](@ref) spends a fixed budget of bond energy to make the target as abundant as possible.

```jldoctest design
julia> using Crafts, Roly, Convex, Clarabel

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false);

julia> ξ, residual = convexdesign(asys, 6; maxε=8.0, optimizer=Clarabel.Optimizer);

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

`maxε` is read according to `εbound`: `:mean` fixes the mean bond energy, `:max` caps the largest one.
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
`εbound` selects whether the mean or the largest bond energy is minimized.

## Designability

[`lineardesign`](@ref) asks a coarser question: is there any choice of parameters at all that isolates the target?
It places the target at zero excess free energy and pushes everything else below.
When that is impossible, the target set is not *designable*, and with `refine=true` the function reports the smallest designable set containing it.

## Polyhedral analysis

The set of parameters that realize a given target is a polyhedral cone, and Crafts can compute its structure exactly with [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html).
These functions work in exact rational arithmetic and need no solver.

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
