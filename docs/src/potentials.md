# Bond potentials

A bond potential says what it costs two bonded particles to move away from the pose they are bonded in.
An entropy solver integrates that cost over all such poses to give a partition function.

## Rigid springs

[`RigidSpringPotential`](@ref) models a bond as a set of harmonic springs between two clouds of patches, one carried by each of the two bonded binding sites.
Only the total stiffness and the spread of the patches matter, so the constructors take just those.

The simplest form spreads patches uniformly along a segment of length `σ` across the site:

```jldoctest pot
julia> using Crafts, Roly

julia> pot = RigidSpringPotential(0.5; k=100.0);
```

Larger `k` or larger `σ` both make the bond stiffer.
In three dimensions the patches cover a `σx` by `σy` sheet:

```jldoctest pot
julia> pot3d = RigidSpringPotential(0.4, 0.6; k=100.0);
```

You can also give the patch positions explicitly, as a `d` by `n` matrix in the binding site's frame:

```jldoctest pot
julia> explicit = RigidSpringPotential([0.0 0.0; -0.25 0.25]; k=100.0);
```

`k` may be a vector of per-patch stiffnesses instead of a number.
To give different binding site colors different patches, pass one matrix per color; see the docstring for details.

## Bond volumes

[`bondvolume`](@ref) is the configuration space volume of a single bond, which is the partition function per bond of a tree:

```jldoctest pot
julia> round(bondvolume(pot); sigdigits=6)
0.109116
```

Two quantities describe how well the two patch clouds fit together.
[`strainenergy`](@ref) is the energy left in the springs when they cannot all relax at once, and [`contactexcess`](@ref) is how far the bond sits above its own minimum when the two sites are placed as they bond.
Both are zero for the constructors above:

```jldoctest pot
julia> strainenergy(pot), contactexcess(pot)
(0.0, 0.0)
```

The Laplace and mean field solvers require `contactexcess` to vanish, and will throw if it does not.

## Solvers

Four solvers turn a potential into a partition function.

| Solver | What it does |
|:---|:---|
| [`TreeLike`](@ref) | Charges `n - 1` bond volumes; exact for structures that really are trees. |
| [`MeanField`](@ref) | Each particle sees its neighbours' patches at their average positions. |
| [`TetheredLaplace`](@ref) | Expands to second order about the bonded configuration, one particle held fixed. |
| [`COMLaplace`](@ref) | The same, in the centre of mass frame. |

On the three-particle chain they differ only for the trimer:

```jldoctest pot
julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> for solver in (TreeLike(), MeanField(), TetheredLaplace(), COMLaplace())
           Ωs = partitionfunctions(AssemblySystem(rules, EntropyModel(pot, solver)))
           println(rpad(nameof(typeof(solver)), 16), round(Ωs[end]; sigdigits=5))
       end
TreeLike        0.07481
MeanField       0.015908
TetheredLaplace 0.07481
COMLaplace      0.07481
```

`TreeLike` also works with no potential at all, in which case every bond contributes its `omega` argument:

```jldoctest pot
julia> model = EntropyModel(TreeLike(0.5));
```

## Custom potentials

Any callable taking two binding site poses works as a potential:

```jldoctest pot
julia> mypotential(pose1, pose2) = 50.0 * sum(abs2, pose2.x - pose1.x);
```

The Laplace solvers will use it directly.
To resolve binding site colors, define [`bondenergy`](@ref) for your type instead.
To make `TreeLike` or `MeanField` work with it you also need [`bondvolume`](@ref), which has no general closed form and so is only defined for `RigidSpringPotential`.

Define [`checkpotential`](@ref) to have a mismatch between a potential and a rule set reported when the [`AssemblySystem`](@ref) is built rather than later.
