# Diffusion and rates

The default rate of `1.0` per reaction sets the timescale arbitrarily.
These functions compute physical rates instead, from how fast the fragments diffuse and how they have to be oriented to bind.

Everything here is in units where the solvent viscosity is one.

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

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()));

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
