abstract type BondPotential end

"""
    RigidSpringPotential

A bond modelled as harmonic springs between two rigid clouds of patches, one carried by each of the
two bonded binding sites,

    V = -ε + Σ_α kα/2 ‖pα₁(q₁) - pα₂(q₂)‖² ,

with the patch `α` of site `i` sitting at `pαᵢ(qᵢ) = xᵢ + R(ψᵢ) aαᵢ`. Call it as
`potential(pose1, pose2, color1, color2)` on the lab poses of the two binding sites; the binding
energy `ε` is not included, as it is carried separately by the reaction network.

Only the total spring constant `k = Σ kα` and the patch (cross-)covariances
`Sᵢⱼ = k⁻¹ Σ kα aαᵢ aαⱼᵀ` enter, so those are what is stored. Since the poses are those of the
binding sites rather than of the particles, and bonded sites coincide, every cloud is taken to be
centred on its own site: patch positions are used only through their spread `aα - ā`.

See [`bondvolume`](@ref) for the configurational integral of a single such bond.
"""
struct RigidSpringPotential{D,F,L} <: BondPotential
    Sij::Matrix{SMatrix{D,D,F,L}}  # S₁₂ of each ordered slot pair; Sij[c2, c1] == Sij[c1, c2]'
    Scc::Vector{SMatrix{D,D,F,L}}  # Sᵢᵢ, the self covariance of each color's own cloud
    k::Matrix{F}
    color_independent::Bool
end

# `k` is the total spring constant; a vector splits it over the patches as the kα.
_patchweights(n, ks::AbstractVector) = (sum(ks), ks ./ sum(ks))
_patchweights(n, k::Real) = (k, fill(inv(n), n))

# With rank(S) < d - 1 the bond keeps a free rotation: Φ_d stays finite, but its stiff limit does
# not, since then some σ̃_μ + σ̃_ν vanishes. That means one patch in 2d, or collinear patches in 3d.
function _checkedcovariance(S::SMatrix{D}) where {D}
    r = rank(S)
    r >= D - 1 || throw(ArgumentError("patch covariance has rank $r, so the bond is free to " *
                                      "rotate; $(D)d needs at least rank $(D - 1)"))
    return S
end

function _compatiblepair(S::SMatrix{D,D,F}, k) where {D,F}
    # Bonded sites face each other (`Roly.standard_offset`), so their relative orientation at
    # contact is R⋆ = standard_rotationᵀ. Compatibility gives S₁₂ = S₁₁ R⋆ and S₂₂ = R⋆ᵀ S₁₁ R⋆,
    # of which only tr(S₂₂) is needed, and that equals tr(S₁₁), so one slot still holds everything.
    Rstar = Roly.standard_rotation(F, Val(D))'
    Sij = fill(SMatrix{D,D,F}(_checkedcovariance(S) * Rstar), 1, 1)
    return RigidSpringPotential{D,F,D * D}(Sij, [S], fill(F(k), 1, 1), true)
end

"""
    RigidSpringPotential(σ=1; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

Two dimensional, with the patches spread uniformly along a segment of length `σ` across the site,
which has second moment `σ²/12`.
"""
RigidSpringPotential(σ::Real=1; k=1) = _compatiblepair(SMatrix{2,2,float(typeof(σ))}(0, 0, 0, σ^2 / 12), k)

"""
    RigidSpringPotential(σx, σy; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

Three dimensional, with the patches spread uniformly over a `σx x σy` sheet across the site, which
has second moments `σx²/12` and `σy²/12`.
"""
RigidSpringPotential(σx::Real, σy::Real; k=1) =
    _compatiblepair(SMatrix{3,3,float(promote_type(typeof(σx), typeof(σy)))}(0, 0, 0, 0,
                                                                            σx^2 / 12, 0, 0, 0,
                                                                            σy^2 / 12), k)

"""
    RigidSpringPotential(A::AbstractMatrix; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

`A` is an explicit `d x n` cloud of patch positions in the binding site's frame, carried by every
site; `k` is either a number or a vector of the per-patch stiffnesses.
"""
function RigidSpringPotential(A::AbstractMatrix; k=1)
    D, n = size(A)
    ktot, w = _patchweights(n, k)
    F = float(promote_type(eltype(A), eltype(w)))
    α = A .- A * w
    return _compatiblepair(SMatrix{D,D,F}(α * (w .* α')), ktot)
end

"""
    RigidSpringPotential(As::AbstractVector{<:AbstractMatrix}; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

One `d x n` cloud of patch positions per binding site color, `As[c]` in the frame of color `c`,
paired column by column; `k` is either a number or a vector of the per-patch stiffnesses. The patches
of two colors are in general not compatible, unless their clouds are related by the rotation that
brings the two sites face to face. If incompatible, the bond carries a [`strainenergy`](@ref).
"""
function RigidSpringPotential(As::AbstractVector{<:AbstractMatrix}; k=1)
    D, n = size(first(As))
    all(A -> size(A) == (D, n), As) ||
        throw(ArgumentError("every color needs a patch cloud of the same size"))
    ktot, w = _patchweights(n, k)
    F = float(promote_type(eltype(w), mapreduce(eltype, promote_type, As)))

    devs = [Matrix{F}(A .- A * w) for A in As]
    nc = length(As)
    Sij = [SMatrix{D,D,F}(devs[c1] * (w .* devs[c2]')) for c1 in 1:nc, c2 in 1:nc]
    foreach(c -> _checkedcovariance(Sij[c, c]), 1:nc)
    return RigidSpringPotential{D,F,D * D}(Sij, [Sij[c, c] for c in 1:nc],
                                           fill(F(ktot), nc, nc), false)
end

# Which storage slots a bond between two colors reads.
_slots(rsp::RigidSpringPotential, color1, color2) = rsp.color_independent ? (1, 1) : (color1, color2)

function (rsp::RigidSpringPotential{D})(p1::Pose{D}, p2::Pose{D}, color1=1, color2=1) where {D}
    c1, c2 = _slots(rsp, color1, color2)
    ΔR = p1.psi' * p2.psi
    return rsp.k[c1, c2] / 2 * (tr(rsp.Scc[c1]) + tr(rsp.Scc[c2]) + normsq(p2.x - p1.x) -
                                2 * tr(rsp.Sij[c1, c2] * ΔR'))
end

"""
    bondenergy(potential, pose1, pose2, color1, color2)

Energy of one bond, given the lab poses and colors of the two binding sites that form it.
"""
bondenergy(potential, pose1, pose2, color1, color2) = potential(pose1, pose2)
bondenergy(rsp::RigidSpringPotential, pose1, pose2, color1, color2) = rsp(pose1, pose2, color1, color2)

"""
    bondvolume(potential, color1=1, color2=1)
    logbondvolume(potential, color1=1, color2=1)

The configurational integral `ω(S, k) = exp(-∑ε) ∫ exp(-V(g)) dg` of a single bond, over the relative pose
`g ∈ SE(d)` of the two binding sites that form it, excluding the binding energy.
"""
bondvolume(potential, color1=1, color2=1) = exp(logbondvolume(potential, color1, color2))

logbondvolume(potential, color1=1, color2=1) = throw(ArgumentError(
    "$(nameof(typeof(potential))) has no closed-form bond volume; use a Laplace solver instead"))

function logbondvolume(rsp::RigidSpringPotential{D}, color1=1, color2=1) where {D}
    c1, c2 = _slots(rsp, color1, color2)
    k = rsp.k[c1, c2]
    σ = _orientedsingularvalues(k * rsp.Sij[c1, c2])
    return D / 2 * log(2π / k) + D * (D - 1) / 4 * log(2π) - _strainenergy(rsp, c1, c2, σ) -
           sum(log(σ[μ] + σ[ν]) for μ in 1:D for ν in (μ + 1):D; init=zero(eltype(σ))) / 2
end

"""
    checkpotential(potential, rules::BindingRules, d::Integer)

Check that `potential` can describe every bond `rules` allows, in `d` dimensions, and throw if it
cannot. 
"""
checkpotential(potential, rules::BindingRules, d::Integer) = nothing

function checkpotential(rsp::RigidSpringPotential{D}, rules::BindingRules, d::Integer) where {D}
    D == d || throw(ArgumentError("$(D)d bond potential cannot describe $(d)d structures"))
    nc = size(rsp.k, 1)
    rsp.color_independent || nc >= ncolors(rules) ||
        throw(ArgumentError("bond potential resolves $nc binding site colors, but the rules use " *
                            "$(ncolors(rules))"))
    return
end

"""
    contactexcess(potential, color1=1, color2=1)

How far above its own minimum a bond sits when its two binding sites are placed as Roly bonds them,
`k [Σ_μ σ̃_μ(S₁₂) - tr(S₁₂ R⋆ᵀ)]`. Zero exactly when the potential is minimised at that pose.

Distinct from [`strainenergy`](@ref), which asks only whether the springs can relax *somewhere*. A
potential is free to relax elsewhere and [`bondvolume`](@ref) will still be right, since it
integrates over all poses; the Laplace solvers expand about the bonded pose and so reject it.
"""
contactexcess(potential, color1=1, color2=1) = 0.0

function contactexcess(rsp::RigidSpringPotential{D,F}, color1=1, color2=1) where {D,F}
    c1, c2 = _slots(rsp, color1, color2)
    S = rsp.Sij[c1, c2]
    # tr(S₁₂ ΔR₁₂ᵀ) is largest over SO(d) at Σ_μ σ̃_μ, so the bonded pose minimises V iff it gets there
    return rsp.k[c1, c2] *
           (sum(_orientedsingularvalues(S)) - tr(S * Roly.standard_rotation(F, Val(D))))
end

"""
    strainenergy(potential, color1=1, color2=1)

Energy `k [tr(S₁₁ + S₂₂)/2 - Σ_μ σ̃_μ(S₁₂)]` left in the springs of a bond that cannot relax them
all at once, in the stiff-spring limit. Non-negative, and zero exactly when the two patch clouds are
compatible.

Note the sign against the binding energies `εs`, which are bond strengths: a bond contributes
`exp(ε - strainenergy)` to a density, so strain has to be *subtracted* from `ε` to fold it in, not
added.
"""
function strainenergy(rsp::RigidSpringPotential, color1=1, color2=1)
    c1, c2 = _slots(rsp, color1, color2)
    return _strainenergy(rsp, c1, c2, _orientedsingularvalues(rsp.k[c1, c2] * rsp.Sij[c1, c2]))
end

# `σ` are the oriented singular values of k*S₁₂, so they already carry the factor of k
_strainenergy(rsp, c1, c2, σ) = rsp.k[c1, c2] / 2 * (tr(rsp.Scc[c1]) + tr(rsp.Scc[c2])) - sum(σ)

# log ω(k, S) for compatible patches of covariance `S`: there the strain cancels and Φ_d leaves
# √(2π/k)^dtot / ∏_{μ<ν} √(λ_μ + λ_ν), with λ the eigenvalues of S.
function _logcompatibleomega(k, S::SMatrix{D,D}) where {D}
    λ = eigvals(Symmetric(S))
    return (D * (D + 1) ÷ 2) / 2 * log(2π / k) -
           sum(log(λ[μ] + λ[ν]) for μ in 1:D for ν in (μ + 1):D) / 2
end

# The oriented singular values: σ̃ = σ except that the smallest one picks up sgn(det A), which is
# what restricts tr(A R) to SO(d) rather than O(d). Everything downstream is symmetric in the σ̃,
# so `svdvals`' descending order is left alone and the last entry is the one that flips.
function _orientedsingularvalues(A::SMatrix{D,D}) where {D}
    σ = svdvals(A)
    return setindex(σ, sign(det(A)) * σ[D], D)
end
