"""
    logdensities(asys::AssemblySystem, ξ)

Compute the log equilibrium number densities of the structures of `asys` as a function of `ξ`, a vector
containing chemical potentials and binding energies.
"""
function logdensities(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return logdensities(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    logdensities(asys::AssemblySystem, ϕs, εs)

Compute the log equilibrium number densities of the structures of `asys` as a function of particle
concentrations (`ϕs`) and binding energies (`εs`).
"""
function logdensities(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return logdensities(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    logyields(asys::AssemblySystem, ξ)

Compute the log equilibrium yields of the structures of `asys` as a function of `ξ`, a vector containing
chemical potentials and binding energies.
"""
function logyields(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return logyields(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    logyields(asys::AssemblySystem, ϕs, εs)

Compute the log equilibrium yields of the structures of `asys` as a function of particle
concentrations (`ϕs`) and binding energies (`εs`).
"""
function logyields(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return logyields(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    logparticledensities(asys::AssemblySystem, ξ)

Compute the log total equilibrium number densities of each particle species used in `asys` as a function
of `ξ`, a vector containing chemical potentials and binding energies.
"""
function logparticledensities(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return logparticledensities(ξ, compositionmatrix(asys), partitionfunctions(asys); ns=nspecies(asys))
end

"""
    logparticledensities(asys::AssemblySystem, ϕs, εs)

Compute the log total equilibrium number densities of each particle species used in `asys` as a
function of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function logparticledensities(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return logparticledensities(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    yields(asys::AssemblySystem, ξ)

Compute the equilibrium yields of the structures of `asys` as a function of `ξ`, a vector containing
chemical potentials and binding energies.
"""
function yields(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return yields(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    yields(asys::AssemblySystem, ϕs, εs)

Compute the equilibrium yields of the structures of `asys` as a function of particle concentrations (`ϕs`)
and binding energies (`εs`).
"""
function yields(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return yields(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    densities(asys::AssemblySystem, ξ)

Compute the equilibrium number densities of the structures of `asys` as a function of `ξ`, a vector containing
chemical potentials and binding energies.
"""
function densities(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return densities(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    densities(asys::AssemblySystem, ϕs, εs)

Compute the equilibrium number densities of the structures of `asys` as a function of particle
concentrations (`ϕs`) and binding energies (`εs`).
"""
function densities(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return densities(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    particledensities(asys::AssemblySystem, ξ)

Compute the total equilibrium number densities of each particle species used in `asys`
as a function of `ξ`, a vector containing chemical potentials and binding energies.
"""
function particledensities(asys::AssemblySystem, ξ)
    return particledensities(ξ, compositionmatrix(asys), partitionfunctions(asys); ns=nspecies(asys))
end

"""
    particledensities(asys::AssemblySystem, ϕs, εs)

Compute the total equilibrium number densities of each particle species used in `asys` as a function
of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function particledensities(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return particledensities(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    chemicalpotentials(asys::AssemblySystem, ϕs, εs)

Compute the chemical potentials of each particle species used in `asys`
as a function of particle concentrations (`ϕs`) and binding energies (`εs`).
"""
function chemicalpotentials(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return chemicalpotentials(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    toyields(densities)

Normalize a list of number densities into yields.
If `densities` has multiple axes, the normalization is carried out over `dims`.
"""
function toyields(densities; dims=1)
    return softmax(log.(abs.(densities)); dims)
end

# The ξ vector realising a given set of particle concentrations and binding energies.
_xivector(ϕs, εs, M, Ωs; solve_kwargs...) = [chemicalpotentials(ϕs, εs, M, Ωs; solve_kwargs...); εs]

function logdensities(ξ, M, Ωs)
    log_ρs = M * ξ .+ log.(Ωs)
    return log_ρs
end
logdensities(ϕs, εs, M, Ωs; kw...) = logdensities(_xivector(ϕs, εs, M, Ωs; kw...), M, Ωs)
densities(ξ, M, Ωs) = exp.(logdensities(ξ, M, Ωs))
densities(ϕs, εs, M, Ωs; kw...) = densities(_xivector(ϕs, εs, M, Ωs; kw...), M, Ωs)

function logyields(ξ, M, Ωs)
    log_ρs = logdensities(ξ, M, Ωs)
    log_ρtot = LogExpFunctions.logsumexp(log_ρs)
    return log_ρs .- log_ρtot
end
logyields(ϕs, εs, M, Ωs; kw...) = logyields(_xivector(ϕs, εs, M, Ωs; kw...), M, Ωs)
yields(ξ, M, Ωs) = exp.(logyields(ξ, M, Ωs))
yields(ϕs, εs, M, Ωs; kw...) = yields(_xivector(ϕs, εs, M, Ωs; kw...), M, Ωs)

function logparticledensities(ξ, M, Ωs; ns=_nspecies(M))
    species = view(M, :, 1:ns)
    log_ρs = logdensities(ξ, M, Ωs)
    return vec(logsumexp(log.(species) .+ log_ρs; dims=1))
end
logparticledensities(ϕs, εs, M, Ωs; kw...) =
    logparticledensities(_xivector(ϕs, εs, M, Ωs; kw...), M, Ωs; ns=length(ϕs))
particledensities(ξ, M, Ωs; ns=_nspecies(M)) = exp.(logparticledensities(ξ, M, Ωs; ns))
particledensities(ϕs, εs, M, Ωs; kw...) = exp.(logparticledensities(ϕs, εs, M, Ωs; kw...))

function chemicalpotentials(ϕs, εs, M, Ωs; alg=nothing, atol=1e-6, rtol=1e-6, maxiters=1_000_000,
                            infval=max(99, 10 * maximum(εs)))
    any(<(0), ϕs) && throw(ArgumentError("Particle concentrations cannot be negative."))

    nμ = length(ϕs)
    nε = length(εs)
    nμ + nε != size(M, 2) && throw(ArgumentError("Lengths of `ϕs` and `εs` does not match the shape of `M`."))

    T = float(promote_type(eltype(ϕs), eltype(εs)))

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
    logNs = [T.(log.(N[mask, k])) for (k, mask) in enumerate(masks)]
    Ms = [hcat(N[mask, :], B[mask, :]) for mask in masks]
    logΩs = [T.(log.(Ωs[mask])) for mask in masks]

    function f!(Δϕ, μs, εs)
        ξ = vcat(μs, εs)
        for (k, i) in enumerate(μrange)
            Δϕ[k] = logsumexp(logNs[k] .+ logΩs[k] .+ Ms[k] * ξ) - log(ϕs[i])
        end
        return Δϕ
    end

    init_μs = fill(T(-1.1 * maximum(εs)), length(μrange))
    prob = NonlinearProblem(f!, init_μs, T.(εs))
    solution = solve(prob, alg; abstol=atol, reltol=rtol, maxiters)

    if solution.retcode == ReturnCode.Stalled
        @warn "Conversion from chemical potentials to particle concentrations stalled. The solution may be inaccurate, proceed with care."
    elseif solution.retcode != ReturnCode.Success
        @error "Conversion from chemical potentials to particle concentrations failed with status $(solution.retcode)."
    end

    μs = fill(T(-infval), nμ)
    μs[μrange] .= solution.u
    return μs
end

function _check_parameterlength(asys::AssemblySystem, ξ)
    n = nspecies(asys) + nbonds(asys)
    length(ξ) != n && throw(ArgumentError("The assembly system takes $n parameters, but `ξ` only has length $(length(ξ))."))
    return
end

function _check_parameterlength(asys::AssemblySystem, ϕs, εs)
    length(ϕs) != nspecies(asys) && throw(ArgumentError("The assembly system contains $(nspecies(asys)) particle species, but `ϕs` only has length $(length(ϕs))."))
    length(εs) != nbonds(asys) && throw(ArgumentError("The assembly system contains $(nbonds(asys)) bond types, but `εs` only has length $(length(εs))."))
    return
end

function _nspecies(M::AbstractMatrix)
    nspc = 0
    for m in eachrow(M)
        if sum(abs, m) == 1
            nspc += 1
        end
    end
    return nspc
end

"""
    density_jacobian(asys::AssemblySystem, ξ)
    density_jacobian(ξ, M, Ωs)

Derivative of the structure number densities with respect to `ξ`, a vector containing chemical
potentials and binding energies.
"""
density_jacobian(ξ, M, Ωs) = densities(ξ, M, Ωs) .* M

function density_jacobian(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return density_jacobian(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    yield_jacobian(asys::AssemblySystem, ξ)
    yield_jacobian(ξ, M, Ωs)

Derivative of the structure yields with respect to `ξ`, a vector containing chemical potentials and
binding energies.
"""
function yield_jacobian(ξ, M, Ωs)
    Ys = yields(ξ, M, Ωs)
    return Ys .* (M .- Ys' * M)
end

function yield_jacobian(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return yield_jacobian(ξ, compositionmatrix(asys), partitionfunctions(asys))
end

"""
    particledensity_jacobian(asys::AssemblySystem, ξ)
    particledensity_jacobian(ξ, M, Ωs; ns=_nspecies(M))

Derivative of the total particle densities with respect to `ξ`, a vector containing chemical
potentials and binding energies.
"""
function particledensity_jacobian(ξ, M, Ωs; ns=_nspecies(M))
    ρN = @view(M[:, 1:ns]) .* densities(ξ, M, Ωs)
    return ρN' * M
end

function particledensity_jacobian(asys::AssemblySystem, ξ)
    _check_parameterlength(asys, ξ)
    return particledensity_jacobian(ξ, compositionmatrix(asys), partitionfunctions(asys); ns=nspecies(asys))
end

function _density_jacobian_ϕ(ρs, M, ns)
    N = @view M[:, 1:ns]
    B = @view M[:, ns+1:end]

    ρN = N .* ρs
    ∂ϕ_∂ξ = ρN' * M
    ∂ϕ_∂μ = @view ∂ϕ_∂ξ[:, 1:ns]
    ∂μ_∂ε = -(∂ϕ_∂μ \ @view ∂ϕ_∂ξ[:, ns+1:end])

    return hcat(ρN / ∂ϕ_∂μ, ρs .* (N * ∂μ_∂ε + B))
end

"""
    density_jacobian(asys::AssemblySystem, ϕs, εs)
    density_jacobian(ϕs, εs, M, Ωs)

Derivative of the structure number densities with respect to particle concentrations (`ϕs`) and
binding energies (`εs`).
"""
function density_jacobian(ϕs, εs, M, Ωs; solve_kwargs...)
    return _density_jacobian_ϕ(densities(ϕs, εs, M, Ωs; solve_kwargs...), M, length(ϕs))
end

function density_jacobian(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return density_jacobian(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end

"""
    yield_jacobian(asys::AssemblySystem, ϕs, εs)
    yield_jacobian(ϕs, εs, M, Ωs)

Derivative of the structure yields with respect to particle concentrations (`ϕs`) and binding
energies (`εs`).
"""
function yield_jacobian(ϕs, εs, M, Ωs; solve_kwargs...)
    log_ρs = logdensities(ϕs, εs, M, Ωs; solve_kwargs...)
    ρs = exp.(log_ρs)
    Ys = exp.(log_ρs .- logsumexp(log_ρs))
    ∂ρs = _density_jacobian_ϕ(ρs, M, length(ϕs))
    return (∂ρs .- Ys .* sum(∂ρs; dims=1)) / sum(ρs)
end

function yield_jacobian(asys::AssemblySystem, ϕs, εs; solve_kwargs...)
    _check_parameterlength(asys, ϕs, εs)
    return yield_jacobian(ϕs, εs, compositionmatrix(asys), partitionfunctions(asys); solve_kwargs...)
end
