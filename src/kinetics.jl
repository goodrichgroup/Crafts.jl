"""
    simulate_kinetics(net::ReactionNetwork, ξ; T, [initial_densities, saveat])

Simulate the assembly kinetics over the reactions of `net` as a function of `ξ`, a vector containing
chemical potentials and binding energies. Unless `initial_densities` is specified, the initial state
consists of monomers at the concentrations implied by `ξ`.

Keyword arguments:
- `T`: total simulation time.
- `initial_densities` (optional): initial densities of all structures at t=0.
- `saveat` (optional): time points at which the simulation state should be stored.
- `sparsejac` (optional): store the Jacobian as a sparse matrix.
- `alg` (optional): the ODE algorithm, `Rodas5P()` by default.

Remaining keyword arguments go to `solve`.

Returns a vector of time points and a matrix of structure number densities of shape `(nstructures, ntimepoints)`.
"""
function simulate_kinetics(net::ReactionNetwork, ξ; T, initial_densities=nothing, saveat=(), kwargs...)
    _check_parameterlength(assemblysystem(net), ξ)
    return _simulate_kinetics(net, ξ; Ts=(0, T), initial_densities, saveat, kwargs...)
end

"""
    simulate_kinetics(net::ReactionNetwork, ϕs, εs; T, kwargs...)

As above, but as a function of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function simulate_kinetics(net::ReactionNetwork, ϕs, εs; solve_kwargs=(;), kwargs...)
    asys = assemblysystem(net)
    _check_parameterlength(asys, ϕs, εs)
    return simulate_kinetics(net, topotentials(asys, ϕs, εs; solve_kwargs...); kwargs...)
end

# One 3x3 block per reaction, on the rows and columns of its `(i, j, k)`: the pattern that
# `stabilitymatrix!` fills. Built from the reactions rather than from a stability matrix, since that
# one needs detailed balance, which the kinetics does not, and since a rate that happens to vanish
# would drop an entry the solver still has to be able to write to.
function _jacobianpattern(rxns, n)
    rows = Vector{Int}(undef, 9 * length(rxns))
    cols = similar(rows)

    at = 0
    for (i, j, k) in rxns, a in (i, j, k), b in (i, j, k)
        rows[at += 1] = a
        cols[at] = b
    end
    return sparse(rows, cols, ones(length(rows)), n, n, (x, y) -> 1.0)
end

_zerojacobian!(J::SparseMatrixCSC) = fill!(nonzeros(J), 0)
_zerojacobian!(J) = fill!(J, 0)

function _make_odefunction_and_ratescale(net::ReactionNetwork, ξ; sparsejac=false)
    rxns = _copy_active_reactions(net)
    kfwd = _copy_active_fwdrates(net)
    kbwd = _copy_active_bwdrates(net)

    log_ρs = logdensities(assemblysystem(net), ξ)
    logsum = zero(eltype(kbwd))

    ndissociating = 0
    for (r, (i, j, k)) in enumerate(rxns)
        log_kbwd = log(kbwd[r]) + log_ρs[i] + log_ρs[j] - log_ρs[k] # stay in log space
        kbwd[r] = log_kbwd
        if isfinite(log_kbwd) # protect against log(0)
            logsum += log_kbwd
            ndissociating += 1
        end
    end

    log_kscale = if ndissociating > 0
        logsum / ndissociating
    elseif !isempty(rxns)
        sum(log, kfwd) / length(kfwd)
    else
        zero(eltype(kfwd))
    end

    # only the backward rates are scaled: with ρ = kscale*u and t = t′/kscale, dρ/dt gains a factor
    # kscale^2, and so does the second-order aggregation term, leaving `kfwd` untouched.
    kbwd .= exp.(kbwd .- log_kscale)
    kscale = exp(log_kscale)
    isfinite(kscale) && kscale > 0 ||
        throw(ArgumentError("rate scale $kscale is out of range; check the kernels and `ξ`"))

    function update_step!(du, u, p, t)
        du .= 0

        for r in eachindex(rxns)
            i, j, k = rxns[r]

            Ka = kfwd[r] * u[i] * u[j]
            Kb = kbwd[r] * u[k]

            Rij = Kb - Ka # We don't need a factor of 2 for i == j, because we add it twice in that case!

            du[i] += Rij
            du[j] += Rij
            du[k] -= Rij
        end
        return
    end

    # Handing the solver the exact Jacobian, rather than letting it difference `update_step!` once per
    # structure, is where nearly all of the time goes otherwise. The derivatives of
    # `Rij = kbwd*u[k] - kfwd*u[i]*u[j]` are the three below; `i == j` needs no special case, since the
    # two contributions are accumulated just as they are in `update_step!`.
    function jacobian!(J, u, p, t)
        _zerojacobian!(J)

        for r in eachindex(rxns)
            i, j, k = rxns[r]

            ∂i = -kfwd[r] * u[j]
            ∂j = -kfwd[r] * u[i]
            ∂k = kbwd[r]

            J[i, i] += ∂i
            J[i, j] += ∂j
            J[i, k] += ∂k

            J[j, i] += ∂i
            J[j, j] += ∂j
            J[j, k] += ∂k

            J[k, i] -= ∂i
            J[k, j] -= ∂j
            J[k, k] -= ∂k
        end
        return
    end

    # The rates are constant in time, so the Rosenbrock methods need not difference for this either.
    function time_gradient!(dT, u, p, t)
        dT .= 0
        return
    end

    jac_prototype = sparsejac ? _jacobianpattern(rxns, nstructures(assemblysystem(net))) : nothing
    return ODEFunction(update_step!; jac=jacobian!, tgrad=time_gradient!, jac_prototype), kscale
end

function _simulate_kinetics(net::ReactionNetwork, ξ; Ts, initial_densities=nothing, alg=Rodas5P(),
                            sparsejac=false, saveat=(), abstol=nothing, solver_kwargs...)
    asys = assemblysystem(net)
    np = nspecies(asys)
    nstr = nstructures(asys)

    f, kscale = _make_odefunction_and_ratescale(net, ξ; sparsejac)
    tscale = inv(kscale)
    ρscale = kscale

    if isnothing(initial_densities)
        initial_densities = vcat(particledensities(asys, ξ), zeros(nstr - np))
    end
    if isnothing(abstol)
        abstol = sum(initial_densities) * 1e-8
    end

    initial_densities = initial_densities / ρscale
    Ts = Ts ./ tscale
    saveat = saveat ./ tscale
    abstol /= ρscale

    prob = ODEProblem(f, initial_densities, Ts)
    sol = solve(prob, alg; saveat, abstol, solver_kwargs...)

    ts = sol.t * tscale
    us = reduce(hcat, sol.u * ρscale)
    return ts, us
end
