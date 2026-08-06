"""
    bondcolors(poly::Polyform)

The binding site colors `(c1, c2)` of every bond of `poly`, in the order [`bonds`](@ref) gives them.
"""
bondcolors(poly::Polyform) = ((_sitecolor(poly, p1, s1), _sitecolor(poly, p2, s2))
                              for ((p1, s1), (p2, s2)) in bonds(poly))

_sitecolor(poly::Polyform, particle, site) =
    Roly.siteloc2color(bindingrules(poly), (poly.particles[particle].species_index, site))

function map_potential(bond_potential, poly::Polyform; embed3d=false)
    sys = bindingrules(poly)
    bs = collect(bonds(poly))

    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    function energy_fn(ξs::AbstractMatrix)
        E = 0.0
        rottype = d == 2 ? Angle2d{eltype(ξs)} : RotXYZ{eltype(ξs)}
        for ((p1, s1), (p2, s2)) in bs
            particle1 = poly.particles[p1]
            particle2 = poly.particles[p2]
            spcs1 = particle1.species_index
            spcs2 = particle2.species_index

            particle_pose1 = embed3d ? Pose{3}(particle1.pose) : particle1.pose
            particle_pose2 = embed3d ? Pose{3}(particle2.pose) : particle2.pose
            site_pose1 = let sp=bindingsites(species(sys, spcs1), s1).pose
                embed3d ? Pose{3}(sp) : sp
            end
            site_pose2 = let sp=bindingsites(species(sys, spcs2), s2).pose
                embed3d ? Pose{3}(sp) : sp
            end

            # transform coordinates such that minimum is at zero.
            pose1 = Pose(ξs[1:d, p1] + particle_pose1.x, rottype(ξs[d+1:dtot, p1]...) * particle_pose1.psi)
            pose2 = Pose(ξs[1:d, p2] + particle_pose2.x, rottype(ξs[d+1:dtot, p2]...) * particle_pose2.psi)

            x1 = pose1 * site_pose1
            x2 = pose2 * site_pose2

            E += bondenergy(bond_potential, x1, x2, _sitecolor(poly, p1, s1),
                            _sitecolor(poly, p2, s2))
        end
        return E
    end

    return energy_fn
end

function hessian(bond_potential, poly; embed3d=false)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    energy_fn = map_potential(bond_potential, poly; embed3d)
    x0 = zeros(numtype(poly), dtot, n)
    return ForwardDiff.hessian(energy_fn, x0 .+ sqrt(eps(eltype(x0))))
end

function orientational_volume(d::Integer)
    d == 2 && return 2π
    d == 3 && return 8π^2
    throw(ArgumentError("orientational volume requires d=2 or d=3"))
end

function inertia_tensor(poly::Polyform; embed3d=false)
    d = dimension(poly)
    com = sum(let x=p.pose.x; x end for p in poly.particles) / nparticles(poly)
    inertia = sum(let x=p.pose.x - com; norm(x)^2 * I(d) - x * x' end for p in poly.particles)
    embed3d || return inertia
    return [inertia zeros(2); zeros(2)' tr(inertia)]
end

abstract type EntropySolver end
struct TetheredLaplace <: EntropySolver end
struct COMLaplace <: EntropySolver end

"""
    MeanFieldSolver()

Entropy solver that is exact at mean field level for a [`RigidSpringPotential`](@ref), where every
particle sees its neighbours' patches at their thermally averaged positions rather than at their
instantaneous ones.

Needs an unstressed structure, so that those averaged positions are the ones the particles are glued
to; each particle then contributes `ω(kᵢ, ⟨Sᵢ⟩)` built from the spread of its own patches, and the
structure's topology enters only through how many patches each particle carries.

Being variational it can only undercount. A dimer it gets exactly right, since tethering one of two
particles leaves a single free one and the product ansatz is then no ansatz at all.
"""
struct MeanFieldSolver <: EntropySolver end

"""
    TreeApproximation(omega=1)

Entropy solver that treats every structure as a tree, charging `n - 1` bond volumes on top of the
structure's orientational volume. A structure that really is a tree has exactly `n - 1` bonds
and the answer is exact. For non-trees, the bond volumes are averaged geometrically over the bonds. 
With a bond potential the bondvolumes are computed from [`bondvolume`](@ref), with no potential each 
bond contributes `omega`.
"""
struct TreeApproximation <: EntropySolver
    omega::Float64
end
TreeApproximation() = TreeApproximation(1.0)

"""
    EntropyModel(potential, solver=COMLaplace(); embed3d=false)
    EntropyModel(solver=TreeApproximation(); embed3d=false)

Everything needed to turn a structure into a partition function: the bond `potential`, the
`solver` used to evaluate it, and whether to embed the structure in 3d. The second form is for
solvers that take no potential.
"""
struct EntropyModel{P,S<:EntropySolver}
    potential::P
    solver::S
    embed3d::Bool
end
EntropyModel(potential, solver::EntropySolver=COMLaplace(); embed3d=false) =
    EntropyModel(potential, solver, embed3d)

function EntropyModel(solver::EntropySolver=TreeApproximation(); embed3d=false)
    solver isa TreeApproximation || throw(ArgumentError(
        "$(nameof(typeof(solver))) needs a bond potential; use `EntropyModel(potential, solver)`."))
    return EntropyModel(nothing, solver, embed3d)
end

"""
    entropy(model::EntropyModel, poly::Polyform)

Partition function of the structure `poly` under `model`.
"""
entropy(m::EntropyModel, poly::Polyform) = _entropy(m.potential, poly, m.solver; m.embed3d)

# A Laplace solver reads off the curvature at the pose Roly bonds the sites in, so it is only an
# expansion of anything if the potential is actually stationary there.
function _checkbondsrelaxed(bond_potential, poly::Polyform)
    for (c1, c2) in bondcolors(poly)
        excess = contactexcess(bond_potential, c1, c2)
        excess > sqrt(eps(typeof(excess))) && throw(ArgumentError(
            "invalid minimizer for Laplace's method: binding site colors ($c1, $c2) bond energy is $excess above the potential's minimum/"))
    end
    return nothing
end

function _entropy(::Nothing, poly::Polyform, solver::TreeApproximation; embed3d)
    d = embed3d ? 3 : dimension(poly)
    return solver.omega^(nparticles(poly) - 1) * orientational_volume(d) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::TreeApproximation; embed3d)
    d = embed3d ? 3 : dimension(poly)
    Ω = orientational_volume(d) / symmetrynumber(poly)
    nparticles(poly) == 1 && return Ω

    # A tree of `n` particles has exactly `n - 1` bonds, so the geometric mean is just their
    # product; anything with cycles is sees `n - 1` "typical" bonds instead.
    logω = mean(logbondvolume(bond_potential, c1, c2) for (c1, c2) in bondcolors(poly))
    return Ω * exp((nparticles(poly) - 1) * logω)
end

_entropy(potential, poly::Polyform, ::MeanFieldSolver; embed3d) = throw(ArgumentError(
    "MeanFieldSolver needs a `RigidSpringPotential`, got $(nameof(typeof(potential)))"))

function _entropy(rsp::RigidSpringPotential, poly::Polyform, ::MeanFieldSolver; embed3d)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    Ω = orientational_volume(d) / symmetrynumber(poly)
    n == 1 && return Ω

    _checkbondsunstressed(rsp, poly)
    # Which particle is tethered changes the answer, so the geometric mean over all n tetherings is
    # taken.
    logω = sum(_logcompatibleomega(k, S) for (k, S) in _patchcovariances(rsp, poly, embed3d))
    return Ω * exp((n - 1) / n * logω)
end

_sitepose(poly::Polyform, particle, site, embed3d) =
    let sp = bindingsites(species(bindingrules(poly), poly.particles[particle].species_index),
                          site).pose
        embed3d ? Pose{3}(sp) : sp
    end

# The total stiffness kᵢ carried by each particle, and Σᵢ, the covariance of every patch on it in
# its own frame. In an unstressed structure a neighbour's averaged patches land exactly on the ones
# facing them, so ⟨Sᵢ⟩ = Σᵢ Rᵢᵀ and Φ_d sees only the eigenvalues of Σᵢ.
function _patchcovariances(rsp::RigidSpringPotential{D,F}, poly::Polyform, embed3d) where {D,F}
    n = nparticles(poly)
    ks = zeros(F, n)
    means = fill(zero(SVector{D,F}), n)
    seconds = fill(zero(SMatrix{D,D,F}), n)

    for ((p1, s1), (p2, s2)) in bonds(poly)
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

function _entropy(bond_potential, poly::Polyform, ::TetheredLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d)
    λs = eigvals(H[dtot+1:end, dtot+1:end])
    S_vib = -0.5 * sum(log, λs / (2π); init=0)
    return orientational_volume(d) * exp(S_vib) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::COMLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    N = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d)
    λs = eigvals(H)
    S_vib = -0.5 * sum(log, λ / (2π) for λ in @view λs[dtot+1:end]; init=0)

    inertia = inertia_tensor(poly; embed3d)
    if d == 2
        detG = N^d * (tr(inertia) + N)
    elseif d == 3
        detG = N^d * det(inertia + N * I(dtot - d))
    else
        throw(ArgumentError("dimension must be 2 or 3"))
    end
    return orientational_volume(d) * sqrt(detG) * exp(S_vib) / symmetrynumber(poly)
end
