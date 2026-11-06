"""
    BondPotential

Supertype of the bond potentials. See [`RigidSpringPotential`](@ref).
"""
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
_patchweights(n, ks::AbstractVector) = (ktot = sum(ks); (ktot, ks ./ ktot))
_patchweights(n, k::Real) = (k, fill(inv(n), n))

# A bond keeps a free rotation when some pair sum λ_μ + λ_ν of the patch covariance vanishes: Φ_d
# stays finite there but its stiff limit does not, and it is exactly those pair sums that the bond
# volume divides by. The smallest of them is λ_{D-1} + λ_D, so the test is on the second smallest
# eigenvalue, measured against the largest to stay scale free. Free rotation means one patch in 2d,
# or collinear patches in 3d.
function _checkedcovariance(S::SMatrix{D}; rtol=sqrt(eps(real(float(eltype(S)))))) where {D}
    λ = eigvals(Symmetric(S))
    λ[2] > rtol * λ[D] || throw(ArgumentError(
        "patch covariance has second smallest eigenvalue $(λ[2]) against largest $(λ[D]), so the " *
        "bond is free to rotate; $(D)d needs at least $(D - 1) stiff directions"))
    return S
end

"""
    RigidSpringPotential{D,F,L}(S, k)

Construct a rigid spring potential with compatible patches, given by the covariance matrix `S`, and
the total stiffness `k`.

Bonded sites face each other (`Roly.standard_offset`), so their relative orientation at contact is
`R⋆ = standard_rotationᵀ`.
"""
function RigidSpringPotential{D,F,L}(S::SMatrix{D,D,F,L}, k::Real) where {D,F,L}
    Rstar = Roly.standard_rotation(F, Val(D))'
    Sij = fill(SMatrix{D,D,F}(_checkedcovariance(S) * Rstar), 1, 1)
    return RigidSpringPotential{D,F,L}(Sij, [S], fill(F(k), 1, 1), true)
end

"""
    RigidSpringPotential(σ=1; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

Two dimensional, with the patches spread uniformly along a segment of length `σ` across the site,
which has second moment `σ²/12`.
"""
function RigidSpringPotential(σ::Real=1; k=1)
    F = float(typeof(σ))
    return RigidSpringPotential{2,F,4}(SMatrix{2,2,F}(0, 0, 0, σ^2 / 12), k)
end

"""
    RigidSpringPotential(σx, σy; k=1)

Harmonic springs of total stiffness `k` between the patch clouds of two bonded binding sites.

Three dimensional, with the patches spread uniformly over a `σx x σy` sheet across the site, which
has second moments `σx²/12` and `σy²/12`.
"""
function RigidSpringPotential(σx::Real, σy::Real; k=1)
    F = float(promote_type(typeof(σx), typeof(σy)))
    return RigidSpringPotential{3,F,9}(SMatrix{3,3,F}(0, 0, 0, 0, σx^2 / 12, 0, 0, 0, σy^2 / 12), k)
end

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
    return RigidSpringPotential{D,F,D * D}(SMatrix{D,D,F}(α * (w .* α')), ktot)
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

A bond potential is called as `potential(pose1, pose2, color1, color2)`, so one that does not resolve
colors still has to accept and ignore them.
"""
bondenergy(potential, pose1, pose2, color1, color2) = potential(pose1, pose2, color1, color2)

# The Hessian of `f` at the bonded configuration `x0`. Bonded sites coincide, so a potential written
# with e.g. `norm(x)` is smooth there as a function but its automatic derivative is 0 * Inf.
function _bondedhessian(f, x0; perturb)
    perturb && return ForwardDiff.hessian(f, x0 .+ sqrt(eps(eltype(x0))))
    H = ForwardDiff.hessian(f, x0)
    all(isfinite, H) || throw(ArgumentError(
        "the Hessian of the bond potential at the bonded configuration is not finite, which happens when " *
        "the potential is not differentiable there as written, e.g. `norm(x)` at `x = 0`. Write it in a " *
        "differentiable form, e.g. `sum(abs2, x)` for `norm(x)^2`, or pass `perturb=true` to " *
        "`EntropyModel` (or `hessian`, `bondvolume`) to evaluate the Hessian a distance √eps away."))
    return H
end

# The two binding site poses of a bond whose second site is displaced by `ξ`, site 1 at the origin
# and site 2 facing it, which is how Roly places bonded sites. `ξ` has the `dtot = d(d+1)/2`
# coordinates of the relative pose the bond volume integrates over, and `ξ = 0` is the configuration
# the Laplace solvers expand about.
function _bondposes(::Val{D}, ξ) where {D}
    F = eltype(ξ)
    rot = D == 2 ? Angle2d{F} : RotXYZ{F}
    p1 = Pose(zero(SVector{D,F}), rot(ntuple(_ -> zero(F), D * (D - 1) ÷ 2)...))
    p2 = Pose(SVector{D,F}(ξ[1:D]), rot(ξ[(D + 1):end]...) * Roly.standard_rotation(F, Val(D)))
    return p1, p2
end

"""
    bondvolume(potential, color1=1, color2=1)
    logbondvolume(potential, color1=1, color2=1)

The configurational integral `ω(S, k) = exp(-∑ε) ∫ exp(-V(g)) dg` of a single bond, over the relative pose
`g ∈ SE(d)` of the two binding sites that form it, excluding the energy at the minimum.

By default this is evaluated using Laplace's approximation. Passing `exact=true` will compute the result exactly, if
the potential implements the exact calculation.

A potential with no closed form is expanded about the bonded pose by automatic differentiation,
which needs the dimension `d` of that pose since only the caller knows it. `perturb=true` evaluates
that Hessian a distance `√eps` away from the bonded pose, for potentials whose automatic derivatives
are not finite there.
"""
bondvolume(potential, color1=1, color2=1; kwargs...) =
    exp(logbondvolume(potential, color1, color2; kwargs...))

function logbondvolume(potential, color1=1, color2=1; d, perturb=false, kwargs...)
    dtot = d * (d + 1) ÷ 2
    V(ξ) = bondenergy(potential, _bondposes(Val(d), ξ)..., color1, color2)
    H = _bondedhessian(V, zeros(dtot); perturb)
    return dtot / 2 * log(2π) - logdet(Symmetric(H)) / 2
end

function logbondvolume(rsp::RigidSpringPotential{D}, color1=1, color2=1;
                       exact=false, d=D, perturb=false) where {D}
    D == d || throw(ArgumentError("$(D)d bond potential asked for a $(d)d bond volume"))
    c1, c2 = _slots(rsp, color1, color2)
    k = rsp.k[c1, c2]
    exact && return D / 2 * log(2π / k) -
                   k * (tr(rsp.Scc[c1]) + tr(rsp.Scc[c2])) / 2 + logphi(k * rsp.Sij[c1, c2])

    σ = _orientedsingularvalues(k * rsp.Sij[c1, c2])
    return D / 2 * log(2π / k) + D * (D - 1) / 4 * log(2π) - _strainenergy(rsp, c1, c2, σ) -
           sum(log(σ[μ] + σ[ν]) for μ in 1:D for ν in (μ + 1):D; init=zero(eltype(σ))) / 2
end

# besselix(ν, x) = besseli(ν, x) exp(-|x|), so this recovers the log of the unscaled function
# without ever forming the overflowing intermediate.
_logbesseli0(x) = log(besselix(0, x)) + abs(x)

"""
    logphi(A)

Return `log Φ_d(A)` with `Φ_d(A) = ∫_SO(d) exp(tr(A R)) dR`, the orientational integral of a rigid
spring bond, exactly.

Closed form in two dimensions, and a one dimensional Bessel integral in three. `Φ_d` grows like
`exp(∑ σ̃_μ)` in the singular values of `A` while the prefactor of ω decays like `exp(-k tr(S₁₁+S₂₂)/2)`,
and the two cancel for compatible patches, so everything is kept in logs.

The three dimensional integrand peaks at the endpoint, ever more sharply as the springs stiffen.
Substituting `t = 1 - u` moves that peak to `t = 0` and pulls its height out as `exp(∑ σ̃_μ)`, leaving
`exp(-(σ̃₁+σ̃₂)t)`, and geometric breakpoints then resolve it however narrow it becomes.
"""
function logphi(A::AbstractMatrix)
    d = size(A, 1)
    if d == 2
        return log(2π) + _logbesseli0(sqrt(tr(A)^2 + (A[1, 2] - A[2, 1])^2))
    elseif d == 3
        σ = _orientedsingularvalues(SMatrix{3,3}(A))
        a = (σ[3] + σ[2]) / 2
        b = (σ[3] - σ[2]) / 2
        s = σ[1] + σ[2]
        g(t) = exp(-s * t) * besselix(0, a * (2 - t)) * besselix(0, b * t)
        edges = unique([0.0; 2 .* exp10.(range(-14, 0; length=29))])
        integral = sum(first(quadgk(g, edges[i], edges[i + 1]; rtol=1e-10))
                       for i in 1:(length(edges) - 1))
        return log(4π^2) + sum(σ) + log(integral)
    end
    throw(ArgumentError("Φ_d is implemented for d = 2 and d = 3, got d = $d"))
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
    contactexcess(rsp::RigidSpringPotential, color1=1, color2=1)

How far above its own minimum a bond sits when its two binding sites are placed as Roly bonds them,
`k [Σ_μ σ̃_μ(S₁₂) - tr(S₁₂ R⋆ᵀ)]`. Zero exactly when the potential is minimised at that pose, so that
`V(contact) = strainenergy + contactexcess`.

Distinct from [`strainenergy`](@ref), which asks only whether the springs can relax *somewhere*. A
potential is free to relax elsewhere and [`bondvolume`](@ref) will still be right, since it
integrates over all poses; the Laplace solvers expand about the bonded pose and so reject it.

Defined for a rigid spring only, where both terms are closed form. Another potential would need its
minimum found numerically, so the Laplace solvers take its expansion point on trust.
"""
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
function _logomega_compatible(k, S::SMatrix{D,D}) where {D}
    λ = eigvals(Symmetric(S))
    return (D * (D + 1) ÷ 2) / 2 * log(2π / k) -
           sum(log(λ[μ] + λ[ν]) for μ in 1:D for ν in (μ + 1):D) / 2
end

# The oriented singular values, ascending: σ̃ = σ except that the smallest picks up sgn(det A),
# which is what restricts tr(A R) to SO(d) rather than O(d). Ascending so `logphi` can index them
# positionally; every other caller uses only symmetric functions of them, so the order is free.
function _orientedsingularvalues(A::SMatrix{D,D}) where {D}
    σ = reverse(svdvals(A))
    return setindex(σ, sign(det(A)) * σ[1], 1)
end

# A Laplace solver reads off the curvature at the pose Roly bonds the sites in, so it is only an
# expansion of anything if the potential is actually stationary there.
# Only a rigid spring reports how far the bonded pose sits above its own minimum, so that is the
# only potential whose expansion point can be checked here rather than taken on trust.
_checkbondsrelaxed(bond_potential, poly::Polyform) = nothing

function _checkbondsrelaxed(rsp::RigidSpringPotential, poly::Polyform)
    for (c1, c2) in bondcolors(poly)
        excess = contactexcess(rsp, c1, c2)
        excess > sqrt(eps(typeof(excess))) && throw(ArgumentError(
            "invalid minimizer for Laplace's method: binding site colors ($c1, $c2) bond energy is $excess above the potential's minimum/"))
    end
    return nothing
end

# Mean field takes the neighbours' patches at their averaged positions, which are the glued ones
# only if every bond really is relaxed where it sits.
function _checkbondsunstressed(rsp::RigidSpringPotential, poly::Polyform)
    _checkbondsrelaxed(rsp, poly)
    for (c1, c2) in bondcolors(poly)
        strain = strainenergy(rsp, c1, c2)
        strain > sqrt(eps(typeof(strain))) && throw(ArgumentError(
            "binding site colors ($c1, $c2) carry a strain energy of $strain, but mean field " *
            "theory is only worked out here for unstressed structures; see `strainenergy`"))
    end
    return nothing
end

_sitepose(poly::Polyform, particle, site, embed3d) =
    let sp = bindingsite(species(bindingrules(poly), poly.particles[particle].speciesindex),
                         site).pose
        embed3d ? Pose{3}(sp) : sp
    end

# The total stiffness kᵢ carried by each particle, and Σᵢ, the covariance of every patch on it in
# its own frame. In an unstressed structure a neighbour's averaged patches land exactly on the ones
# facing them, so ⟨Sᵢ⟩ = Σᵢ Rᵢᵀ and Φ_d sees only the eigenvalues of Σᵢ.
function _patchcovariances(rsp::RigidSpringPotential{D,F}, poly::Polyform, embed3d;
                           bondlist=bonds(poly)) where {D,F}
    n = nparticles(poly)
    ks = zeros(F, n)
    means = fill(zero(SVector{D,F}), n)
    seconds = fill(zero(SMatrix{D,D,F}), n)

    for ((p1, s1), (p2, s2)) in bondlist
        slots = _slots(rsp, _sitecolor(poly, p1, s1), _sitecolor(poly, p2, s2))
        k = rsp.k[slots...]
        for (p, s, slot) in ((p1, s1, slots[1]), (p2, s2, slots[2]))
            pose = _sitepose(poly, p, s, embed3d)
            ks[p] += k
            means[p] += k * pose.x
            # the cloud is centred on its own site, so the site sits at its mean patch position and
            # carrying it over to the particle's frame is the parallel axis theorem
            seconds[p] += k * (pose.psi * rsp.Scc[slot] * pose.psi' + pose.x * pose.x')
        end
    end
    return ((ks[i], seconds[i] / ks[i] - (means[i] / ks[i]) * (means[i] / ks[i])') for i in 1:n)
end

