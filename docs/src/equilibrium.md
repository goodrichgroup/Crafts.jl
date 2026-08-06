# Equilibrium

At equilibrium the number density of a structure is set by its composition, its partition function, and a vector `ξ` holding one chemical potential per species and one binding energy per bond type.

## From concentrations to potentials

Most of the time you know the total particle concentrations `ϕs` and the bond energies `εs`, not the chemical potentials.
[`topotentials`](@ref) solves for the chemical potentials that reproduce those concentrations and returns the full `ξ`:

```jldoctest eq
julia> using Crafts, Roly

julia> rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare);

julia> asys = AssemblySystem(rules, EntropyModel(TreeLike()));

julia> ϕs, εs = fill(0.1, 3), fill(8.0, 2);

julia> ξ = topotentials(asys, ϕs, εs);

julia> length(ξ) == nspecies(asys) + nbonds(asys)
true
```

[`chemicalpotentials`](@ref) returns just the chemical potential part.

Everything below takes either `ξ` or the pair `ϕs, εs`, and the two forms are interchangeable.

## Densities and yields

[`densities`](@ref) gives the number density of each structure, [`yields`](@ref) the fraction of structures of each kind, and [`particledensities`](@ref) the total density of each species across all structures:

```jldoctest eq
julia> round.(densities(asys, ϕs, εs); digits=4)
6-element Vector{Float64}:
 0.0135
 0.0018
 0.0135
 0.0117
 0.0117
 0.0748

julia> round.(yields(asys, ϕs, εs); digits=4)
6-element Vector{Float64}:
 0.1063
 0.0144
 0.1063
 0.092
 0.092
 0.5891

julia> round.(particledensities(asys, ϕs, εs); digits=4)
3-element Vector{Float64}:
 0.1
 0.1
 0.1
```

The particle densities come back as the `ϕs` that went in, which is what `topotentials` solved for.

Each has a `log` counterpart — [`logdensities`](@ref), [`logyields`](@ref), [`logparticledensities`](@ref) — which is the better choice when the densities span many orders of magnitude.

## Derivatives

The Jacobians give the response of the densities to the parameters:

```jldoctest eq
julia> size(yield_jacobian(asys, ϕs, εs))
(6, 5)
```

Rows are structures.
When called with `ξ`, columns are the entries of `ξ`; when called with `ϕs, εs`, columns are the concentrations followed by the bond energies.
[`density_jacobian`](@ref) and [`particledensity_jacobian`](@ref) work the same way.
