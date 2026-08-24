# Diffusion and rates

The default rate of `1.0` per reaction sets the timescale arbitrarily.
These functions compute physical rates instead, from how fast the fragments diffuse and how they have to be oriented to bind.

Everything here is in units where the solvent viscosity is one.

## Dimension

Every rate here is three-dimensional.
The mobilities are Stokes' law `1/(6πr)` and the Rotne-Prager-Yamakawa tensor built from the 3d Oseen kernel, the friction tensor carries the six rigid-body motions of a 3d body, and [`smoluchowskirate`](@ref) is the steady-state flux onto an absorbing sphere in three dimensions.
None of these have a two-dimensional counterpart.
A disk in an unbounded two-dimensional Stokes flow has no finite mobility at all, since the flow cannot be brought to rest at infinity (Stokes' paradox), and the absorbing-disk problem has no steady state, so its capture rate keeps decaying like `4πD/log(4Dt/R²)` rather than settling on a rate constant.

A two-dimensional structure is therefore handled as a flat rigid cluster of spheres suspended in a three-dimensional fluid: its coordinates are padded with `z = 0`, and the rate comes out the same whether or not the entropy model embeds.
The partition functions do not follow along, so an [`EntropyModel`](@ref) built without `embed3d=true` scores those structures in 2d while these kernels move them in 3d.
The two meet in the dissociation rates, which [`simulate_kinetics`](@ref) and [`stabilitymatrix`](@ref) get from the forward rates by detailed balance, so the mismatch is not a matter of an overall constant.
[`ReactionNetwork`](@ref) warns when [`entropydimension`](@ref) is not 3.

## Hydrodynamics

A structure is treated as a rigid cluster of spheres.
[`frictiontensor`](@ref) assembles the `6x6` grand friction tensor from the positions and radii, using Rotne-Prager-Yamakawa pair mobilities:

```jldoctest rates
julia> using Crafts, Roly

julia> θ = frictiontensor([[0.0, 0, 0], [1.2, 0, 0]], [0.5, 0.5]);

julia> size(θ)
(6, 6)
```

Blocks are translation-translation, translation-rotation, and rotation-rotation, about the origin.

[`diffusiontensor`](@ref) inverts it about the cluster's centre of diffusion, and [`diffusionconstants`](@ref) returns the orientationally averaged translational and rotational diffusion constants:

```jldoctest rates
julia> round.(diffusionconstants([[0.0, 0, 0], [1.2, 0, 0]], [0.5, 0.5]); sigdigits=5)
(0.075157, 0.092968)
```

A `Polyform` can be passed directly, in which case the sphere radii are taken from how far each species' binding sites sit from its centre:

```jldoctest rates
julia> rules = BindingRules([1 1 2 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike(); embed3d=true));

julia> round.(diffusionconstants(structures(asys)[end]); sigdigits=5)
(0.079577, 0.098524)
```

## Binding rates

Two kernels are provided, both usable directly with [`rate!`](@ref).

[`smoluchowskirate`](@ref) is the diffusion-limited rate for two spheres that react on contact whatever their orientation, `4πDR`:

```jldoctest rates
julia> net = ReactionNetwork(asys; fwdkernel=smoluchowskirate);

julia> round.(fwdrates(net); sigdigits=5)
1-element Vector{Float64}:
 2.6667
```

[`orientedbindingrate`](@ref) also requires the two partners to be correctly oriented, with reactive patches of a given size:

```jldoctest rates
julia> net = ReactionNetwork(asys; fwdkernel=orientedbindingrate);

julia> round.(fwdrates(net); sigdigits=5)
1-element Vector{Float64}:
 0.31601
```

The orientational requirement costs about a factor of eight here.
`siteradius` sets how large the reactive patches are, and a larger one always gives a faster rate:

```jldoctest rates
julia> kernel(rxn) = orientedbindingrate(rxn; siteradius=1.0);

julia> net = ReactionNetwork(asys; fwdkernel=kernel);

julia> round.(fwdrates(net); sigdigits=5)
1-element Vector{Float64}:
 0.73293
```

Widening the patches until they cover the spheres recovers the Smoluchowski rate.
Both kernels are also callable with explicit parameters rather than a reaction; see their docstrings.
