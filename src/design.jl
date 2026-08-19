"""
    lineardesign(M, idxs; optimizer, kwargs...)
    lineardesign(asys::AssemblySystem, idxs; optimizer, kwargs...)

Find parameters that make the target structures `idxs` the only ones at zero excess free energy.

  - `idxs`: index or indices of the target structures
  - `optimizer`: any `MathOptInterface` optimizer, e.g. `Clarabel.Optimizer`. Configure solver
    tolerances on the optimizer itself, with `MOI.OptimizerWithAttributes`
  - `preprocess = true`: drop structures built from particles or bonds the targets do not use
  - `refine = true`: if the targets are not designable on their own, optimize for the smallest
    designable set containing them
  - `atol = 1e-6`: tolerance for deciding that a structure sits at zero excess free energy
  - `silent = true`: suppress solver output
  - `infval = 100`: finite stand-in for infinite parameters

Returns `(ξ, residual)`. Requires `Convex` to be loaded.
"""
function lineardesign end

"""
    maxyielddesign(M, idxs; omegas, maxdensity, optimizer, kwargs...)
    maxyielddesign(asys::AssemblySystem, idxs; maxdensity, optimizer, kwargs...)

Maximize the yield of the target structures `idxs` within a bond energy budget.

  - `idxs`: index or indices of the target structures
  - `maxdensity`: cap on the total particle density, in the units `omegas` is expressed in
  - `relative_yields = nothing`: densities to hold the targets at, relative to each other, one entry
    per target; `nothing` holds them equal
  - `omegas`: entropic partition function of every structure
  - `energy_budget = 1`: cap on the bond energies, applied to `energy_measure` of them
  - `energy_measure = mean`: what the budget bounds, e.g. `mean` or `maximum`
  - `uniform_energy = false`: hold every bond type at the same energy
  - `target_stoichiometry = false`: keep the particle densities in the targets' stoichiometric ratio
  - `optimizer`: any `MathOptInterface` optimizer, e.g. `Clarabel.Optimizer`
  - `preprocess`, `silent`, `infval`: as in [`lineardesign`](@ref)

The budget is a cap, not a target: bond energy also favours the larger off-target structures, so the
solution can come in under it.

`target_stoichiometry` is a cap too, since holding the ratios exactly is not a convex constraint. It
caps each species at the targets' ratio, which leaves under-supply allowed; that is expected to bind
rather than to relax, and a warning reports the solves where it does not.

Returns `(ξ, residual)`, where `residual` is `log Σ_{j∉idxs} ρ_j/ρ_i`. Requires `Convex` to be loaded.
"""
function maxyielddesign end

"""
    minenergydesign(M, idxs; minyield, omegas, maxdensity, optimizer, kwargs...)
    minenergydesign(asys::AssemblySystem, idxs; minyield, maxdensity, optimizer, kwargs...)

Find the weakest bonds that still reach a yield of `minyield` for the target structures `idxs`.

The inverse of [`maxyielddesign`](@ref): the yield enters as a constraint and the bond energies
become the objective.

  - `minyield`: yield the targets must reach, as a fraction of all structures present
  - `energy_measure = mean`: what is minimized, e.g. `mean` or `maximum` of the bond energies
  - `uniform_energy = false`: hold every bond type at the same energy
  - `target_stoichiometry = false`: as in [`maxyielddesign`](@ref), and a cap there too
  - other arguments as in [`maxyielddesign`](@ref), `energy_budget` excepted: here the energy is the
    objective, not a constraint

Returns `(ξ, residual)`, where `residual` is the achieved `energy_measure`. Requires `Convex` to be
loaded.
"""
function minenergydesign end

for f in (:lineardesign, :maxyielddesign, :minenergydesign)
    @eval function $f(args...; kwargs...)
        throw(ArgumentError("`" * $(string(f)) * "` needs Convex.jl and a solver. Run " *
                            "`using Convex, Clarabel` and pass `optimizer=Clarabel.Optimizer`."))
    end
end

lineardesign(asys::AssemblySystem, idxs; kwargs...) =
    lineardesign(compositionmatrix(asys), idxs; kwargs...)
maxyielddesign(asys::AssemblySystem, idxs; kwargs...) =
    maxyielddesign(compositionmatrix(asys), idxs; omegas=partitionfunctions(asys), kwargs...)
minenergydesign(asys::AssemblySystem, idxs; kwargs...) =
    minenergydesign(compositionmatrix(asys), idxs; omegas=partitionfunctions(asys), kwargs...)

# Drops structures that contain bonds or particles absent from the targets, since no choice of
# parameters can suppress them relative to the targets.
function _preprocessdesign(M, idxs)
    nstructs, npars = size(M)

    D = M[idxs, :]
    element_mask = vec(sum(D; dims=1) .!= 0)
    all(element_mask) && return ones(Bool, nstructs), ones(Bool, npars), idxs

    structure_mask = vec((M * .!element_mask) .== 0)
    new_idxs = [sum(structure_mask[1:i]) for i in idxs]

    return structure_mask, element_mask, new_idxs
end
