module CraftsConvexExt

using Crafts, Convex, Statistics
using Crafts: _nspecies, _preprocessdesign, infapprox

_asindices(idxs) = idxs isa Number ? [Int(idxs)] : collect(Int, idxs)

# A failed solve still leaves values on the variable, and they generally violate the constraints, so
# the status has to be checked before those values are handed back.
function _checkstatus(problem, what)
    status = problem.status
    status == Convex.MOI.OPTIMAL && return nothing
    if status in (Convex.MOI.INFEASIBLE, Convex.MOI.INFEASIBLE_OR_UNBOUNDED)
        throw(ArgumentError("no parameters satisfy the $what problem; the requested yields or " *
                            "concentration bounds cannot be met at once"))
    elseif status == Convex.MOI.DUAL_INFEASIBLE
        throw(ArgumentError("the $what problem is unbounded; constrain the bond energies or the " *
                            "total concentration"))
    end
    @warn "The $what solver returned $status; the result may be inaccurate, proceed with care."
    return nothing
end

# Everything maxyielddesign and minenergydesign share: the preprocessed system, the ratio objective
# ρ_j/ρ_i for non-targets, the total particle density, and the relative-yield equalities.
function _designterms(M, idxs, omegas, relative_yields, preprocess)
    idxs = _asindices(idxs)
    nμ = _nspecies(M)

    if preprocess
        structure_mask, element_mask, idxs = _preprocessdesign(M, idxs)
        M = M[structure_mask, element_mask]
        omegas = omegas[structure_mask]
        missing_pars = findall(.!element_mask)
        nμ -= count(<=(nμ), missing_pars)
    else
        missing_pars = Int[]
    end

    i = first(idxs)
    others = setdiff(axes(M, 1), idxs)

    A = M[others, :] .- M[i, :]'
    s = omegas[others] ./ omegas[i]
    ns = vec(sum(M[:, 1:nμ]; dims=2))

    rest = idxs[2:end]
    dMs = M[[i], :] .- M[rest, :]
    rs = log.(relative_yields[1] .* omegas[rest]) .- log.(relative_yields[2:end] .* omegas[i])

    return (; M, omegas, nμ, idxs, A, s, ns, dMs, rs, missing_pars)
end

# The bond energies as a vector of scalar expressions rather than the slice `x[(nμ + 1):end]`.
# Convex cannot iterate an expression, so `mean` of a slice fails while `mean` of this vector builds
# the same affine atom; `maximum` and anything else the caller passes work on either.
_bondenergies(x, nμ, npars) = [x[i] for i in (nμ + 1):npars]

# Ties every bond type to a single energy by equating neighbouring entries, which stays affine. With
# one bond type there is nothing to tie.
function _adduniform!(constraints, x, nμ, npars)
    npars > nμ + 1 && push!(constraints, x[(nμ + 2):npars] == x[(nμ + 1):(npars - 1)])
    return constraints
end

function _restore(ξ, missing_pars, infval)
    for m in missing_pars
        insert!(ξ, m, -Inf)
    end
    return infapprox(ξ, infval)
end

function Crafts.lineardesign(M::AbstractMatrix, idxs; optimizer, preprocess=true, refine=true,
                             atol=1e-6, silent=true, infval=100)
    idxs = _asindices(idxs)
    nμ = _nspecies(M)

    if preprocess
        structure_mask, element_mask, idxs = _preprocessdesign(M, idxs)
        M = M[structure_mask, element_mask]
        missing_pars = findall(.!element_mask)
        nμ -= count(<=(nμ), missing_pars)
    else
        missing_pars = Int[]
    end

    npars = size(M, 2)
    if npars > 1
        x = Variable(npars)
        A = M[setdiff(axes(M, 1), idxs), :]
        B = M[idxs, :]
        problem = minimize(maximum(A * x), B * x == 0,
                           sum(x[(nμ + 1):end]) - sum(x[1:nμ]) == npars)
        Convex.solve!(problem, optimizer; silent)
        _checkstatus(problem, "lineardesign")
        ξ = vec(x.value)
        residual = problem.optval
    else
        ξ = [0.0]
        residual = -Inf
    end

    # Structures that also sit at zero excess free energy cannot be suppressed; if the caller allows
    # it, design for the smallest designable set that contains the targets instead.
    chimeras = setdiff(findall(isapprox.(M * ξ, 0.0; atol=atol)), idxs)
    if refine && !isempty(chimeras)
        @warn "Target set is not designable, optimizing for minimal enclosing designable set..."
        ξ, residual = Crafts.lineardesign(M, union(idxs, chimeras); optimizer, preprocess=false,
                                          refine=false, atol, silent, infval)
    end

    return _restore(ξ, missing_pars, infval), residual
end

function Crafts.maxyielddesign(M::AbstractMatrix, idxs; omegas, relative_yields=nothing,
                               energy_budget=1, maxdensity=1, energy_measure=mean,
                               uniform_energy=false, optimizer, preprocess=true, silent=true,
                               infval=100)
    nidx = idxs isa Number ? 1 : length(idxs)
    relative_yields = something(relative_yields, ones(nidx))
    t = _designterms(M, idxs, omegas, relative_yields, preprocess)

    npars = size(t.M, 2)
    if npars > 1
        x = Variable(npars)
        constraints = Convex.Constraint[Convex.logsumexp(t.M * x + log.(t.omegas) + log.(t.ns)) <= log(maxdensity),
                                        energy_measure(_bondenergies(x, t.nμ, npars)) <= energy_budget]
        uniform_energy && _adduniform!(constraints, x, t.nμ, npars)
        isempty(t.rs) || push!(constraints, t.dMs * x == t.rs)

        problem = minimize(Convex.logsumexp(t.A * x + log.(t.s)), constraints)
        Convex.solve!(problem, optimizer; silent)
        _checkstatus(problem, "maxyielddesign")
        ξ = vec(x.value)
        residual = problem.optval
    else
        ξ = [0.0]
        residual = -Inf
    end

    return _restore(ξ, t.missing_pars, infval), residual
end

function Crafts.minenergydesign(M::AbstractMatrix, idxs; minyield, omegas, relative_yields=nothing,
                                maxdensity=1, energy_measure=mean, uniform_energy=false, optimizer,
                                preprocess=true, silent=true, infval=100)
    nidx = idxs isa Number ? 1 : length(idxs)
    relative_yields = something(relative_yields, ones(nidx))

    # The relative-yield equalities pin ρ_j/ρ_i, so the targets hold a fixed multiple K of ρ_i and
    # their combined yield is K/(K+R). A floor on that yield is an upper bound on R.
    K = sum(relative_yields) / relative_yields[1]
    0 < minyield < 1 || throw(ArgumentError("`minyield` must lie in (0, 1)."))

    t = _designterms(M, idxs, omegas, relative_yields, preprocess)

    npars = size(t.M, 2)
    if npars > 1
        x = Variable(npars)
        constraints = Convex.Constraint[Convex.logsumexp(t.M * x + log.(t.omegas) + log.(t.ns)) <= log(maxdensity),
                                        Convex.logsumexp(t.A * x + log.(t.s)) <=
                                        log(K * (1 - minyield) / minyield)]
        uniform_energy && _adduniform!(constraints, x, t.nμ, npars)
        isempty(t.rs) || push!(constraints, t.dMs * x == t.rs)

        problem = minimize(energy_measure(_bondenergies(x, t.nμ, npars)), constraints)
        Convex.solve!(problem, optimizer; silent)
        _checkstatus(problem, "minenergydesign")
        ξ = vec(x.value)
        residual = problem.optval
    else
        ξ = [0.0]
        residual = -Inf
    end

    return _restore(ξ, t.missing_pars, infval), residual
end

end # module CraftsConvexExt
