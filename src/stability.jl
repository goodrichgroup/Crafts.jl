"""
    StabilityBasis

The basis a stability matrix is expressed in. [`DensityBasis`](@ref) is the plain one, in structure
number densities; [`SymmetricBasis`](@ref) is the similarity transform that makes it symmetric.
"""
abstract type StabilityBasis end

"""
    DensityBasis()

Coordinates are structure number densities, so `S = ∂(dρ/dt)/∂ρ` as it stands.
"""
struct DensityBasis <: StabilityBasis end

"""
    SymmetricBasis()

Coordinates are rescaled by `D = diagm(sqrt.(ρeq))`, giving `inv(D) * S * D`. Detailed balance makes
that symmetric, so its spectrum is real, and it is the same spectrum as in [`DensityBasis`](@ref)
since the two differ by a similarity transform.
"""
struct SymmetricBasis <: StabilityBasis end

# Mixing the Float64 rates with whatever `ξ` is, which is how ForwardDiff gets through.
_stabilitytype(net::ReactionNetwork, ξ) = promote_type(eltype(ξ), eltype(fwdrates(net)))

"""
    stabilitymatrix(net::ReactionNetwork, ξ, [basis]; scale=1)
    stabilitymatrix!(S, net::ReactionNetwork, ξ, [basis]; scale=1)

The `nstructures x nstructures` stability matrix of `net` about the equilibrium at `ξ`, expressed in
`basis` (default [`DensityBasis`](@ref)). `scale` multiplies every rate, and so the whole matrix.

Requires [`isdetailedbalanced`](@ref): the equilibrium densities are only a steady state of the kinetics
when each reaction's forward and backward rates agree, so there is nothing to linearise about otherwise.
"""
function stabilitymatrix(net::ReactionNetwork, ξ, basis::StabilityBasis=DensityBasis(); kwargs...)
    n = nstructures(net)
    return stabilitymatrix!(zeros(_stabilitytype(net, ξ), n, n), net, ξ, basis; kwargs...)
end

stabilitymatrix!(S, net::ReactionNetwork, ξ; kwargs...) =
    stabilitymatrix!(S, net, ξ, DensityBasis(); kwargs...)

"""
    stabilityjacobian(net::ReactionNetwork, ξ, [basis]; scale=1)
    stabilityjacobian!(J, net::ReactionNetwork, ξ, [basis]; scale=1)

Derivative of [`stabilitymatrix`](@ref) with respect to `ξ`, of shape
`nstructures x nstructures x length(ξ)`. Requires [`isdetailedbalanced`](@ref), as above.
"""
function stabilityjacobian(net::ReactionNetwork, ξ, basis::StabilityBasis=DensityBasis(); kwargs...)
    n = nstructures(net)
    return stabilityjacobian!(zeros(_stabilitytype(net, ξ), n, n, length(ξ)), net, ξ, basis; kwargs...)
end

stabilityjacobian!(J, net::ReactionNetwork, ξ; kwargs...) =
    stabilityjacobian!(J, net, ξ, DensityBasis(); kwargs...)

# Inactive reactions are skipped in the loop rather than filtered out up front, so that nothing here
# allocates per call and nothing goes stale against a later `rate!`.
function stabilitymatrix!(S, net::ReactionNetwork, ξ, ::DensityBasis; scale=1)
    asys = assemblysystem(net)
    _check_parameterlength(asys, ξ)
    isdetailedbalanced(net) ||
        throw(ArgumentError("stability matrix needs a detailed balanced network; see `isdetailedbalanced`"))
    ρeq = densities(ξ, compositionmatrix(asys), partitionfunctions(asys))

    S .= 0
    rxns, kfwd, kbwd = reactions(net), fwdrates(net), bwdrates(net)
    for r in eachindex(rxns)
        net.active[r] || continue
        i, j, k = rxns[r]

        Ka = kfwd[r] * scale
        Kb = (ρeq[i] * ρeq[j] / ρeq[k]) * kbwd[r] * scale

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

function stabilityjacobian!(J, net::ReactionNetwork, ξ, ::DensityBasis; scale=1)
    asys = assemblysystem(net)
    _check_parameterlength(asys, ξ)
    isdetailedbalanced(net) ||
        throw(ArgumentError("stability jacobian needs a detailed balanced network; see `isdetailedbalanced`"))
    M = compositionmatrix(asys)
    ρeq = densities(ξ, M, partitionfunctions(asys))
    ∂ρeq = permutedims(ρeq .* M)   # `density_jacobian` would evaluate the densities a second time

    J .= 0
    scratch = zeros(_stabilitytype(net, ξ), length(ξ))

    rxns, kfwd, kbwd = reactions(net), fwdrates(net), bwdrates(net)
    for r in eachindex(rxns)
        net.active[r] || continue
        i, j, k = rxns[r]

        Ka = kfwd[r] * scale
        Kb = kbwd[r] * scale

        @views @. begin
            scratch = -Ka * ∂ρeq[:, j]
            J[i, i, :] += scratch
            J[j, i, :] += scratch
            J[k, i, :] -= scratch

            scratch = -Ka * ∂ρeq[:, i]
            J[i, j, :] += scratch
            J[j, j, :] += scratch
            J[k, j, :] -= scratch

            scratch = Kb * (-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] +
                            ∂ρeq[:, i] * ρeq[j] / ρeq[k] +
                            ∂ρeq[:, j] * ρeq[i] / ρeq[k])
            J[i, k, :] += scratch
            J[j, k, :] += scratch
            J[k, k, :] -= scratch
        end
    end
    return J
end

function stabilitymatrix!(S, net::ReactionNetwork, ξ, ::SymmetricBasis; scale=1)
    asys = assemblysystem(net)
    _check_parameterlength(asys, ξ)
    isdetailedbalanced(net) ||
        throw(ArgumentError("stability matrix needs a detailed balanced network; see `isdetailedbalanced`"))
    ρeq = densities(ξ, compositionmatrix(asys), partitionfunctions(asys))
    ρeq_sqrt = sqrt.(ρeq)

    S .= 0
    rxns, kfwd, kbwd = reactions(net), fwdrates(net), bwdrates(net)
    for r in eachindex(rxns)
        net.active[r] || continue
        i, j, k = rxns[r]

        Ka = kfwd[r] * scale
        Kb = kbwd[r] * scale

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

function stabilityjacobian!(J, net::ReactionNetwork, ξ, ::SymmetricBasis; scale=1)
    asys = assemblysystem(net)
    _check_parameterlength(asys, ξ)
    isdetailedbalanced(net) || throw(ArgumentError("stability jacobian needs a detailed balanced network; see `isdetailedbalanced`"))
    M = compositionmatrix(asys)
    ρeq = densities(ξ, M, partitionfunctions(asys))
    ρeq_sqrt = sqrt.(ρeq)
    ∂ρeq = permutedims(ρeq .* M)   # `density_jacobian` would evaluate the densities a second time

    J .= 0
    scratch = zeros(_stabilitytype(net, ξ), length(ξ))

    rxns, kfwd, kbwd = reactions(net), fwdrates(net), bwdrates(net)
    for r in eachindex(rxns)
        net.active[r] || continue
        i, j, k = rxns[r]

        Ka = kfwd[r] * scale
        Kb = kbwd[r] * scale

        @views @. begin
            scratch = -Ka * ((ρeq_sqrt[j] / ρeq_sqrt[i]) * ∂ρeq[:, i] +
                             (ρeq_sqrt[i] / ρeq_sqrt[j]) * ∂ρeq[:, j]) / 2
            J[i, j, :] += scratch
            J[j, i, :] += scratch

            # ik
            scratch = (ρeq_sqrt[i] / ρeq_sqrt[k] * ∂ρeq[:, j] +
                       ρeq[j] / (2ρeq_sqrt[i] * ρeq_sqrt[k]) * ∂ρeq[:, i] -
                       ρeq_sqrt[i] * ρeq[j] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])
            J[i, k, :] += Kb * scratch
            J[k, i, :] += Ka * scratch

            # jk
            scratch = (ρeq_sqrt[j] / ρeq_sqrt[k] * ∂ρeq[:, i] +
                       ρeq[i] / (2ρeq_sqrt[j] * ρeq_sqrt[k]) * ∂ρeq[:, j] -
                       ρeq_sqrt[j] * ρeq[i] / (2ρeq_sqrt[k]^3) * ∂ρeq[:, k])
            J[j, k, :] += Kb * scratch
            J[k, j, :] += Ka * scratch

            J[i, i, :] += -Ka * ∂ρeq[:, j]
            J[j, j, :] += -Ka * ∂ρeq[:, i]
            J[k, k, :] += -Kb * (-ρeq[i] * ρeq[j] / ρeq[k]^2 * ∂ρeq[:, k] +
                                 ∂ρeq[:, i] * ρeq[j] / ρeq[k] +
                                 ∂ρeq[:, j] * ρeq[i] / ρeq[k])
        end
    end
    return J
end

"""
    correlationtime(net::ReactionNetwork, ξ; nconserved=nspecies(assemblysystem(net)), scale=1)
    correlationtime(S::AbstractMatrix; nconserved)

The slowest relaxation time of the linearised kinetics, `-1/λ` for the least negative eigenvalue of the
stability matrix that is not one of its zero modes.

- `nconserved`: how many zero eigenvalues to skip, one per conserved particle count.
- `rtol`: eigenvalues within `rtol * maximum(abs, S)` of zero count as zero modes.

Returns `Inf` when the selected mode is itself a zero mode, meaning the kinetics has more conserved
quantities than `nconserved` says.

The network method builds `S` in the [`SymmetricBasis`](@ref), which makes the spectrum real; that
spectrum is shared with [`DensityBasis`](@ref), so the basis only affects accuracy. It inherits the
[`isdetailedbalanced`](@ref) requirement from [`stabilitymatrix`](@ref).
"""
function correlationtime(net::ReactionNetwork, ξ; nconserved=nspecies(assemblysystem(net)),
                         rtol=1e-12, kwargs...)
    return correlationtime(Symmetric(stabilitymatrix(net, ξ, SymmetricBasis(); kwargs...));
                           nconserved, rtol)
end

function correlationtime(S::Symmetric; nconserved, rtol=1e-12)
    C = maximum(abs, S)
    iszero(C) && return Inf
    return _slowestmode(eigvals(Symmetric(S ./ C)) .* C, nconserved, rtol * C)
end

function correlationtime(S::AbstractMatrix; nconserved, rtol=1e-12)
    C = maximum(abs, S)
    iszero(C) && return Inf
    return _slowestmode(sort!(real.(schur(S ./ C).values) .* C), nconserved, rtol * C)
end

function _slowestmode(λs, nconserved, atol)
    n = length(λs)
    n > nconserved ||
        throw(ArgumentError("$n modes cannot include $nconserved conserved quantities and still relax"))
    λ = λs[n-nconserved]
    return λ < -atol ? -inv(λ) : Inf
end
