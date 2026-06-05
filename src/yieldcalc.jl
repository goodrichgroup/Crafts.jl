"""
    yields(strs::StructureCollection, ξ)

Compute the equilibrium yields of the structures in the collection `strs` as a function of `ξ`, a vector containing 
chemical potentials and binding energies.
"""
function yields(strs::StructureCollection, ξ)
    _check_parameterlength(strs, ξ)
    return yields(ξ, strs.M, strs.Zs)
end

"""
    yields(strs::StructureCollection, ϕs, εs)

Compute the equilibrium yields of the structures in the collection `strs` as a function of particle concentrations (`ϕs`)
and binding energies (`εs`).
"""
function yields(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return yields(ϕs, εs, strs.M, strs.Zs)
end

"""
    densities(strs::StructureCollection, ξ)

Compute the equilibrium number densities of the structures in the collection `strs` as a function of `ξ`, a vector containing 
chemical potentials and binding energies.
"""
function densities(strs::StructureCollection, ξ)
    _check_parameterlength(strs, ξ)
    return densities(ξ, strs.M, strs.Zs)
end

"""
    densities(strs::StructureCollection, ϕs, εs)

Compute the equilibrium number densities of the structures in the collection `strs` as a function of particle 
concentrations (`ϕs`) and binding energies (`εs`).
"""
function densities(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return densities(ϕs, εs, strs.M, strs.Zs)
end

"""
    particle_densities(strs::StructureCollection, ξ)

Compute the total equilibrium number densities of each particle species used in the structure collection `strs`
as a function of `ξ`, a vector containing chemical potentials and binding energies.
"""
function particle_densities(strs::StructureCollection, ξ)
    return particle_densities(ξ, strs.M, strs.Zs)
end

"""
    chemical_potentials(strs::StructureCollection, ϕs, εs)

Compute the chemical potentials of each particle species of the used in the structure collection `strs`
as a function of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function chemical_potentials(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return chemical_potentials(ϕs, εs, strs.M, strs.Zs)
end

"""
    toyields(densities)

Normalize a list of number densities into yields.
If `densities` has multiple axes, the normalization is carried out over `dims`.
"""
function toyields(densities; dims=1)
    return softmax(log.(abs.(densities)); dims)
end


function logdensities(ξ, M, Zs)
    log_ρs = M * ξ .+ log.(Zs)
    return log_ρs
end
densities(ξ, M, Zs) = exp.(logdensities(ξ, M, Zs))
densities(ϕs, εs, M, Zs; solve_kwargs...) = densities([chemical_potentials(ϕs, εs, M, Zs; solve_kwargs...); εs], M, Zs)

function _particle_densities(ξ, M, ns, Zs)
    return ns' * densities(ξ, M, Zs)
end
particle_densities(ξ, M, Zs) = _particle_densities(ξ, M, view(M, :, 1:n_species(M)), Zs)
particle_densities(ϕs, εs, M, Zs; solve_kwargs...) = particle_densities([chemical_potentials(ϕs, εs, M, Zs; solve_kwargs...); εs], M, Zs)

function chemical_potentials(ϕs, εs, M, Zs; atol=1e-6, rtol=1e-6, maxiters=1_000_000)
    if any(<(-atol), ϕs)
        throw(ArgumentError("Particle concentrations cannot be negative."))
    end
    nμ = length(ϕs)
    N = M[:, 1:nμ]
    B = M[:, nμ+1:end]

    logNs = [log.(N[N[:, i] .> 0, i]) for i in 1:nμ]
    Ns = [N[N[:, i] .> 0, :] for i in 1:nμ]
    Bs = [B[N[:, i] .> 0, :] for i in 1:nμ]
    logZs = [log.(Zs[N[:, i] .> 0]) for i in 1:nμ]

    logϕplus1(μs, εs, i) = LogExpFunctions.logsumexp([0; logNs[i] .+ logZs[i] .+ Ns[i]*μs .+ Bs[i]*εs])
    f(μs, εs) = [logϕplus1(μs, εs, i) - log(ϕs[i] + 1) for i in 1:nμ] 
    
    init_μs = -1.1 * mean(εs) * ones(nμ)
    prob = NonlinearProblem(f, init_μs, εs, abstol=atol, reltol=rtol)
    solution = solve(prob, TrustRegion(); maxiters)

    if solution.retcode == ReturnCode.Success
        return Vector(solution.u)
    elseif solution.retcode == ReturnCode.Stalled
        @warn "solution status stalled, proceed with care"
        return Vector(solution.u)
    else
        @error "solution status $(solution.retcode)"
        return fill(Missing, nμ)
    end
end

function logyields(ξ, M, Zs)
    log_ρs = logdensities(ξ, M, Zs)
    log_ρtot = LogExpFunctions.logsumexp(log_ρs)
    return log_ρs .- log_ρtot
end
yields(ξ, M, Zs) = exp.(logyields(ξ, M, Zs))
yields(ϕs, εs, M, Zs; solve_kwargs...) = yields([chemical_potentials(ϕs, εs, M, Zs; solve_kwargs...); εs], M, Zs)


function ∂ρ∂μ(ξ, M, Zs)
    ρs = densities(ξ, M, Zs)
    return ρs .* M
end
function ∂Y∂μ(ξ, M, Zs)
    ρs = densities(ξ, M, Zs)
    Ys = yields(ξ, M, Zs) # just to be safe from numerical issues, recompute
    Σρ = sum(ρs)

    ∂ρs = ∂ρ∂μ(ξ, M, Zs)
    return (∂ρs - Ys .* sum(∂ρs, dims=1)) / Σρ
end
function ∂ρ∂ϕ(ϕs, εs, M, Zs)
    np = length(ϕs)
    N = @view M[:, 1:np]
    B = @view M[:, np+1:end]

    ρs = densities(ϕs, εs, M, Zs)
    ∂ϕ∂μ =  N' * Diagonal(ρs) * N
    ∂ϕ∂ε =  N' * Diagonal(ρs) * B

    ∂μ∂ϕ = inv(∂ϕ∂μ)
    ∂μ∂ε = -∂μ∂ϕ * ∂ϕ∂ε

    ∂ρ∂ϕ = ρs .* N * ∂μ∂ϕ
    ∂ρ∂ε = ρs .* (N * ∂μ∂ε + B)

    return hcat(∂ρ∂ϕ, ∂ρ∂ε)
end

function ∂Y∂ϕ(ϕs, εs, M, Zs)
    ρs = densities(ϕs, εs, M, Zs)
    Ys = yields(ϕs, εs, M, Zs) # just to be safe from numerical issues, recompute
    Σρ = sum(ρs)

    ∂ρs = ∂ρ∂ϕ(ϕs, εs, M, Zs)
    return (∂ρs - Ys .* sum(∂ρs, dims=1)) / Σρ
end

function ∂ϕ∂μ(ξ, M, Zs)
    np = n_species(M)
    N = @view M[:, 1:np]
    B = @view M[:, np+1:end]

    ρs = densities(ξ, M, Zs)
    ∂ϕ∂μ =  N' * Diagonal(ρs) * N
    ∂ϕ∂ε =  N' * Diagonal(ρs) * B
    return hcat(∂ϕ∂μ, ∂ϕ∂ε)
end

function _check_parameterlength(strs, ξ)
    length(ξ) != size(strs.M, 2) && throw(ArgumentError("The assembly system takes $(size(strs.M, 2)) parameters, but `ξ` only has length $(length(ξ))."))
    return
end

function _check_parameterlength(strs, ϕs, εs)
    length(ϕs) != size(strs.sys)[1] && throw(ArgumentError("The assembly system contains $(size(strs.sys)[1]) particle species, but `ϕs` only has length $(length(ϕs))."))
    length(εs) != size(strs.sys)[2] && throw(ArgumentError("The assembly system contains $(size(strs.sys)[2]) bond types, but `εs` only has length $(length(εs))."))
    return
end