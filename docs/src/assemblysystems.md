# Assembly systems

An [`AssemblySystem`](@ref) is a set of binding rules together with a way of scoring each structure those rules allow.

## Binding rules

Rules come from Roly.
Each row of the matrix says that a site on one species binds a site on another, as `[species1 site1 species2 site2]`:

```jldoctest asys
julia> using Crafts, Roly

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
2d BindingRules[n=3, k=2]
```

`UnitSquare` is one of Roly's particle shapes; its four sites are numbered anticlockwise.
See the Roly documentation for the other shapes and for how to build your own.

## Entropy models

An [`EntropyModel`](@ref) pairs a bond potential with a solver that turns it into a partition function per structure.
The simplest model ignores the potential entirely and gives every bond the same weight:

```jldoctest asys
julia> model = EntropyModel(TreeLike());
```

With a potential, pass it first:

```jldoctest asys
julia> model_stiff = EntropyModel(RigidSpringPotential(0.5; k=100.0), TetheredLaplace());
```

See [Bond potentials](potentials.md) for the potentials and the four solvers.
Pass `embed3d=true` to treat a two-dimensional structure as a rigid body in three dimensions.

## Building the system

```jldoctest asys
julia> asys = AssemblySystem(rules, model)
AssemblySystem[n=3, k=2]
```

Nothing is enumerated until you ask for it.
The structures, the composition matrix and the partition functions are computed on first use and cached:

```jldoctest asys
julia> nstructures(asys)
6

julia> compositionmatrix(asys)
6×5 Matrix{Int64}:
 1  0  0  0  0
 0  1  0  0  0
 0  0  1  0  0
 1  1  0  1  0
 0  1  1  0  1
 1  1  1  1  1

julia> round.(partitionfunctions(asys); digits=4)
6-element Vector{Float64}:
 6.2832
 6.2832
 6.2832
 6.2832
 6.2832
 6.2832
```

Each row of the composition matrix is one structure.
The first `nspecies` columns count particles of each species, the rest count bonds of each type.
Structures are ordered by size.

[`structures`](@ref) returns the Roly `Polyform`s themselves, which carry the geometry:

```jldoctest asys
julia> structures(asys)[end]
Polyform{2}[n=3]

julia> nbonds(structures(asys)[end])
2
```

## Limiting the enumeration

Many rule sets allow far more structures than you want, or infinitely many.
`maxsize` caps the number of particles per structure:

```jldoctest asys
julia> small = AssemblySystem(rules, model; maxsize=2, verbose=false);

julia> nstructures(small)
5

julia> iscomplete(small)
false
```

A capped system warns the first time it enumerates, since everything computed from it is a truncation; `verbose=false` silences that, as above.
[`iscomplete`](@ref) reports whether the enumeration finished, and is `missing` before anything has been enumerated.

To find out how many structures a rule set allows without paying for the enumeration, use [`countstructures`](@ref):

```jldoctest asys
julia> countstructures(asys)
PolyformCount[n=6, exact, largest=3]
```

It returns an exact count when the search is cheap and a statistical estimate otherwise.
