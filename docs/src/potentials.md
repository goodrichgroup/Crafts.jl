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

## Harmonic bonds

[`HarmonicPotential`](@ref) models a bond by a quadratic form `V = qᵀ H q / 2` in the relative pose `q = (Δx, φ)` of the two binding sites, about the pose they are bonded in.
`Δx` is the position of the second site and `φ` the rotation vector that carries it from its facing orientation, both in the frame of the first site.
To first order `φ` equals the `RotXYZ` angles the Laplace solvers displace a site by, so `H` is exactly the Hessian they see.
In three dimensions `H` is a `6 x 6` matrix with the translations first:

```jldoctest pot
julia> using LinearAlgebra

julia> harmonic = HarmonicPotential(diagm([100.0, 100.0, 100.0, 5.0, 2.0, 8.0]));
```

The translations are those of the site origins, which coincide at contact.
A Hessian measured for another reference point on the second particle differs by a shear of its translation-rotation blocks and has to be carried over first.
`H` also has to be unchanged when the two sites are exchanged, since Roly may list either one first; `symmetrize=true` projects it onto that subspace.

The Hessian of any other potential at contact gives the harmonic bond that all Laplace solvers treat identically to it:

```jldoctest pot
julia> fromsprings = HarmonicPotential(pot3d; d=3);

julia> bondvolume(fromsprings) ≈ bondvolume(pot3d)
true
```

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
`TreeLike` and `MeanField` expand it about the bonded pose by automatic differentiation, for which [`bondvolume`](@ref) needs the dimension `d` of the bond.
[`HarmonicPotential`](@ref), given the potential and `d`, stores that expansion once.

Define [`checkpotential`](@ref) to have a mismatch between a potential and a rule set reported when the [`AssemblySystem`](@ref) is built rather than later.
