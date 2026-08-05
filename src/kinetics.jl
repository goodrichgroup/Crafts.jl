"""
    simulate_kinetics(asys::AssemblySystem, ξ; T, kernel, [initial_densities=nothing, saveat=[]])

Simulate the assembly kinetics of the structures of `asys` as a function of
`ξ`, a vector containing chemical potentials and binding energies. Unless `initial_densities` is specified,
the initial state consists of monomers at the concentrations implied by `ξ`.

Keyword arguments:
- `T`: total simulation time.
- `kernel`: aggregation kernel used to build the reaction network.
- `initial_densities` (optional): initial densities of all structures at t=0.
- `saveat` (optional): time points at which the simulation state should be stored.

Returns a vector of time points and a matrix of structure number densities of shape `(nstructures, ntimepoints)`.
"""
function simulate_kinetics(asys::AssemblySystem, ξ; T, initial_densities=nothing, saveat=[], kernel, kwargs...)
    _check_parameterlength(asys, ξ)
    ts, ρs = _simulate_kinetics(asys, ξ; Ts=(0,T), kernel, initial_densities, saveat, kwargs...)
    return ts, ρs
end
function simulate_kinetics(asys::AssemblySystem, ϕs, εs; kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return simulate_kinetics(asys, [chemicalpotentials(asys, ϕs, εs); εs]; kwargs...)
end


function kinetic_network(asys::AssemblySystem, ξ; maxbonds, kernel, brkkernel=kernel)
    reactions, ks, fs = generate_reactionnetwork(structures(asys); maxlevel=maxbonds, aggkernel=kernel, brkkernel)

    log_ρs = logdensities(asys, ξ)
    fs = [fs[r] * exp(log_ρs[i] + log_ρs[j] - log_ρs[k]) for (r, (i,j,k)) in enumerate(reactions)]
    kscale = exp(mean(log.(fs))) # geometrical mean
    fs /= kscale

    function update_step!(du, u, p, t)
        du .= 0

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * u[i] * u[j]
            Kb = fs[r] * u[k]

            Rij = Kb - Ka # We don't need a factor of 2 for i == j, because we add it twice in that case!

            du[i] += Rij
            du[j] += Rij
            du[k] -= Rij
        end
        return
    end
    return update_step!, kscale
end

function _simulate_kinetics(asys::AssemblySystem, ξ; Ts, kernel, brkkernel=kernel, maxbonds=Inf, initial_densities=nothing, alg=Rodas5(), saveat=[], abstol=nothing, solver_kwargs...)
    np = nspecies(asys)
    nstr = nstructures(asys)

    step, kscale = kinetic_network(asys, ξ; kernel, brkkernel, maxbonds)
    tscale = inv(kscale)
    ρscale = kscale

    if isnothing(initial_densities)
        initial_densities = vcat(particledensities(asys, ξ), zeros(nstr - np))
    end
    if isnothing(abstol)
        abstol = sum(initial_densities) * 1e-8
    end
    
    initial_densities /= ρscale
    Ts = Ts ./ tscale
    saveat = saveat ./ tscale
    abstol /= ρscale

    prob = ODEProblem(step, initial_densities, Ts)
    sol = solve(prob, alg; saveat, abstol, solver_kwargs...)

    ts = sol.t * tscale
    us = reduce(hcat, sol.u * ρscale)
    return ts, us
end

function stability_matrix(asys::AssemblySystem; symmetrize=false, kernel, brkkernel=kernel, maxbonds)
    reactions, ks, fs = generate_reactionnetwork(structures(asys); maxlevel=maxbonds, aggkernel=kernel, brkkernel)
    M = compositionmatrix(asys)
    Ωs = partitionfunctions(asys)

    function Sfn!(S, ξ)
        S .= 0
        ρeq = densities(ξ, M, Ωs)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r]
            Kb = (ρeq[i] * ρeq[j] / ρeq[k]) * fs[r]

            S[i, i] += -Ka * ρeq[j] 
            S[i, j] += -Ka * ρeq[i] 
            S[i, k] += Kb

            S[j, i] += -Ka * ρeq[j] 
            S[j, j] += -Ka * ρeq[i] 
            S[j, k] += Kb

            S[k, i] += Ka * ρeq[j] 
            S[k, j] += Ka * ρeq[i] 
            S[k, k] += -Kb
        end
        return S 
    end

    function ∂S∂μfn!(S, ξ)
        S .= 0
        ρeq = densities(ξ, M, Ωs)
        ∂ρeq = permutedims(∂ρ∂μ(ξ, M, Ωs))

        scratch = zero(ξ)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r]
            Kb = fs[r] 

            @views @. begin 
                scratch = -Ka * ∂ρeq[:, j] 
                S[i, i, :] += scratch
                S[j, i, :] += scratch
                S[k, i, :] -= scratch


                scratch = -Ka * ∂ρeq[:, i] 
                S[i, j, :] += scratch
                S[j, j, :] += scratch
                S[k, j, :] -= scratch

                scratch = Kb *(-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] + 
                            ∂ρeq[:, i] * ρeq[j] / ρeq[k] 
                            + ∂ρeq[:, j] * ρeq[i] / ρeq[k])
                S[i, k, :] += scratch
                S[j, k, :] += scratch
                S[k, k, :] -= scratch
            end
        end
        return S
    end

    #####################
    #### Return inv(D) S D, where D = sqrt(ρᵢ)
    function Ssym_fn!(S, ξ; scale=1)
        ρeq = densities(ξ, M, Ωs)
        ρeq_sqrt = sqrt.(ρeq)

        S .= 0
        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * scale
            Kb = fs[r] * scale

            ij = ρeq_sqrt[i] * ρeq_sqrt[j] 
            ik = ρeq_sqrt[i] / ρeq_sqrt[k] * ρeq[j]
            jk = ρeq_sqrt[j] / ρeq_sqrt[k] * ρeq[i]

            S[i, i] += -Ka * ρeq[j]
            S[i, j] += -Ka * ij
            S[i, k] += Kb * ik

            S[j, i] += -Ka * ij
            S[j, j] += -Ka * ρeq[i]
            S[j, k] += Kb * jk

            S[k, i] += Ka * ik
            S[k, j] += Ka * jk
            S[k, k] += -Kb * (ρeq[i] * ρeq[j] / ρeq[k])
        end
        return S
    end

    function ∂Ssym∂μ_fn!(S, ξ; scale=1)
        ρeq = densities(ξ, M, Ωs)
        ρeq_sqrt = sqrt.(ρeq)
        ∂ρeq = permutedims(∂ρ∂μ(ξ, M, Ωs))

        S .= 0
        scratch = zero(ξ)

        for r in eachindex(reactions)
            i, j, k = reactions[r]

            Ka = ks[r] * scale
            Kb = fs[r] * scale

            @views @. begin
                scratch = -Ka * ((ρeq_sqrt[j] / ρeq_sqrt[i]) * ∂ρeq[:, i] + (ρeq_sqrt[i] / ρeq_sqrt[j]) * ∂ρeq[:, j]) / 2
                S[i, j, :] += scratch
                S[j, i, :] += scratch

                # ik
                scratch = (ρeq_sqrt[i] / ρeq_sqrt[k] * ∂ρeq[:, j] + ρeq[j] / (2ρeq_sqrt[i] * ρeq_sqrt[k]) * ∂ρeq[:, i] -
                    ρeq_sqrt[i] * ρeq[j] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])
                S[i, k, :] += Kb * scratch
                S[k, i, :] += Ka * scratch

                # jk
                scratch = (ρeq_sqrt[j] / ρeq_sqrt[k] * ∂ρeq[:, i] + ρeq[i] / (2ρeq_sqrt[j] * ρeq_sqrt[k]) * ∂ρeq[:, j] -
                ρeq_sqrt[j] * ρeq[i] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])

                S[j, k, :] += Kb * scratch
                S[k, j, :] += Ka * scratch

                S[i, i, :] += -Ka * ∂ρeq[:, j] 
                S[j, j, :] += -Ka * ∂ρeq[:, i] 
                S[k, k, :] += -Kb * (-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] + ∂ρeq[:, i] * ρeq[j] / ρeq[k] + ∂ρeq[:, j] * ρeq[i] / ρeq[k])
            end
        end
        return S
    end

    if symmetrize
        return Ssym_fn!, ∂Ssym∂μ_fn!
    else
        return Sfn!, ∂S∂μfn!
    end
end

function τc(S; np)
    ns = size(S, 1)

    C = maximum(abs, S)
    S = S / C

    λ = partialsort!(schur(S).values * C, ns-np, by=real)
    return -inv(real(λ))
end
