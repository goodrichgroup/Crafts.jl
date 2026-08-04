"""
    yields(strs::StructureCollection, ξ)

Compute the equilibrium yields of the structures in the collection `strs` as a function of `ξ`, a vector containing 
chemical potentials and binding energies.
"""
function yields(strs::StructureCollection, ξ)
    _check_parameterlength(strs, ξ)
    return yields(ξ, strs.M, strs.Ωs)
end

"""
    yields(strs::StructureCollection, ϕs, εs)

Compute the equilibrium yields of the structures in the collection `strs` as a function of particle concentrations (`ϕs`)
and binding energies (`εs`).
"""
function yields(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return yields(ϕs, εs, strs.M, strs.Ωs)
end

"""
    densities(strs::StructureCollection, ξ)

Compute the equilibrium number densities of the structures in the collection `strs` as a function of `ξ`, a vector containing 
chemical potentials and binding energies.
"""
function densities(strs::StructureCollection, ξ)
    _check_parameterlength(strs, ξ)
    return densities(ξ, strs.M, strs.Ωs)
end

"""
    densities(strs::StructureCollection, ϕs, εs)

Compute the equilibrium number densities of the structures in the collection `strs` as a function of particle 
concentrations (`ϕs`) and binding energies (`εs`).
"""
function densities(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return densities(ϕs, εs, strs.M, strs.Ωs)
end

"""
    particledensities(strs::StructureCollection, ξ)

Compute the total equilibrium number densities of each particle species used in the structure collection `strs`
as a function of `ξ`, a vector containing chemical potentials and binding energies.
"""
function particledensities(strs::StructureCollection, ξ)
    return particledensities(ξ, strs.M, strs.Ωs)
end

"""
    chemicalpotentials(strs::StructureCollection, ϕs, εs)

Compute the chemical potentials of each particle species of the used in the structure collection `strs`
as a function of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function chemicalpotentials(strs::StructureCollection, ϕs, εs)
    _check_parameterlength(strs, ϕs, εs)
    return chemicalpotentials(ϕs, εs, strs.M, strs.Ωs)
end

"""
    toyields(densities)

Normalize a list of number densities into yields.
If `densities` has multiple axes, the normalization is carried out over `dims`.
"""
function toyields(densities; dims=1)
    return softmax(log.(abs.(densities)); dims)
end

function logdensities(ξ, M, Ωs)
    log_ρs = M * ξ .+ log.(Ωs)
    return log_ρs
end
densities(ξ, M, Ωs) = exp.(logdensities(ξ, M, Ωs))
densities(ϕs, εs, M, Ωs; solve_kwargs...) = densities([chemicalpotentials(ϕs, εs, M, Ωs; solve_kwargs...); εs], M, Ωs)

function logyields(ξ, M, Ωs)
    log_ρs = logdensities(ξ, M, Ωs)
    log_ρtot = LogExpFunctions.logsumexp(log_ρs)
    return log_ρs .- log_ρtot
end
yields(ξ, M, Ωs) = exp.(logyields(ξ, M, Ωs))
yields(ϕs, εs, M, Ωs; solve_kwargs...) = yields([chemicalpotentials(ϕs, εs, M, Ωs; solve_kwargs...); εs], M, Ωs)

function logparticledensities(ξ, M, Ωs)
    ns = view(M, :, 1:nspecies(M))
    log_ρs = logdensities(ξ, M, Ωs)
    return vec(logsumexp(log.(ns) .+ log_ρs; dims=1))
end
particledensities(ξ, M, Ωs) = exp.(logparticledensities(ξ, M, Ωs))
particledensities(ϕs, εs, M, Ωs; solve_kwargs...) = particledensities([chemicalpotentials(ϕs, εs, M, Ωs; solve_kwargs...); εs], M, Ωs)

function chemicalpotentials(ϕs, εs, M, Ωs; atol=1e-6, rtol=1e-6, maxiters=1_000_000, infval=max(99, 10 * maximum(εs)))
    any(<(0), ϕs) && throw(ArgumentError("Particle concentrations cannot be negative."))

    nμ = length(ϕs)
    nε = length(εs)
    nμ + nε != size(M, 2) && throw(ArgumentError("Lengths of `ϕs` and `εs` does not match the shape of `M`."))

    μrange = findall(>(0), ϕs)
    if length(μrange) != nμ
        μforbid = setdiff(1:nμ, μrange)
        strrange = findall(vec(all(iszero, M[:, μforbid]; dims=2)))
    else
        strrange = axes(M, 1)
    end

    N = M[strrange, μrange]
    B = M[strrange, nμ+1:end]
    Ωs = Ωs[strrange]

    masks = [N[:, k] .> 0 for k in eachindex(μrange)]
    logNs = [log.(N[mask, k]) for (k, mask) in enumerate(masks)]
    Ms = [hcat(N[mask, :], B[mask, :]) for mask in masks]
    logΩs = [log.(Ωs[mask]) for mask in masks]

    function f!(Δϕ, μs, εs)
        ξ = vcat(μs, εs)
        for (k, i) in enumerate(μrange)
            Δϕ[k] = logsumexp(logNs[k] .+ logΩs[k] .+ Ms[k] * ξ) - log(ϕs[i])
        end
        return Δϕ
    end

    init_μs = -1.1 * maximum(εs) * ones(length(μrange))
    prob = NonlinearProblem(f!, init_μs, εs, abstol=atol, reltol=rtol)
    solution = solve(prob; maxiters)

    if solution.retcode == ReturnCode.Stalled
        @warn "Conversion from chemical potentials to particle concentrations stalled. The solution may be inaccurate, proceed with care."
    elseif solution.retcode != ReturnCode.Success
        @error "Conversion from chemical potentials to particle concentrations failed with status $(solution.retcode)."
    end

    μs = fill(-infval, nμ)
    μs[μrange] .= solution.u
    return μs
end

function _check_parameterlength(strs, ξ)
    length(ξ) != size(strs.M, 2) && throw(ArgumentError("The assembly system takes $(size(strs.M, 2)) parameters, but `ξ` only has length $(length(ξ))."))
    return
end

function _check_parameterlength(strs, ϕs, εs)
    length(ϕs) != nspecies(strs.sys) && throw(ArgumentError("The assembly system contains $(nspecies(strs.sys)) particle species, but `ϕs` only has length $(length(ϕs))."))
    length(εs) != nbonds(strs.sys) && throw(ArgumentError("The assembly system contains $(nbonds(strs.sys)) bond types, but `εs` only has length $(length(εs))."))
    return
end

function nspecies(M::AbstractMatrix)
    nspc = 0
    for m in eachrow(M)
        if sum(abs, m) == 1
            nspc += 1
        end
    end
    return nspc
end





# function ∂ρ∂μ(ξ, M, Ωs)
#     ρs = densities(ξ, M, Ωs)
#     return ρs .* M
# end
# function ∂Y∂μ(ξ, M, Ωs)
#     ρs = densities(ξ, M, Ωs)
#     Ys = yields(ξ, M, Ωs) # just to be safe from numerical issues, recompute
#     Σρ = sum(ρs)

#     ∂ρs = ∂ρ∂μ(ξ, M, Ωs)
#     return (∂ρs - Ys .* sum(∂ρs, dims=1)) / Σρ
# end
# function ∂ρ∂ϕ(ϕs, εs, M, Ωs)
#     np = length(ϕs)
#     N = @view M[:, 1:np]
#     B = @view M[:, np+1:end]

#     ρs = densities(ϕs, εs, M, Ωs)
#     ∂ϕ∂μ =  N' * Diagonal(ρs) * N
#     ∂ϕ∂ε =  N' * Diagonal(ρs) * B

#     ∂μ∂ϕ = inv(∂ϕ∂μ)
#     ∂μ∂ε = -∂μ∂ϕ * ∂ϕ∂ε

#     ∂ρ∂ϕ = ρs .* N * ∂μ∂ϕ
#     ∂ρ∂ε = ρs .* (N * ∂μ∂ε + B)

#     return hcat(∂ρ∂ϕ, ∂ρ∂ε)
# end

# function ∂Y∂ϕ(ϕs, εs, M, Ωs)
#     ρs = densities(ϕs, εs, M, Ωs)
#     Ys = yields(ϕs, εs, M, Ωs) # just to be safe from numerical issues, recompute
#     Σρ = sum(ρs)

#     ∂ρs = ∂ρ∂ϕ(ϕs, εs, M, Ωs)
#     return (∂ρs - Ys .* sum(∂ρs, dims=1)) / Σρ
# end

# function ∂ϕ∂μ(ξ, M, Ωs)
#     np = nspecies(M)
#     N = @view M[:, 1:np]
#     B = @view M[:, np+1:end]

#     ρs = densities(ξ, M, Ωs)
#     ∂ϕ∂μ =  N' * Diagonal(ρs) * N
#     ∂ϕ∂ε =  N' * Diagonal(ρs) * B
#     return hcat(∂ϕ∂μ, ∂ϕ∂ε)
# end
