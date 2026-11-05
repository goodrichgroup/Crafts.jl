"""
    bondcolors(poly::Polyform)

The binding site colors `(c1, c2)` of every bond of `poly`, in the order `Roly.bonds` gives them.
"""
bondcolors(poly::Polyform) = ((_sitecolor(poly, p1, s1), _sitecolor(poly, p2, s2))
                              for ((p1, s1), (p2, s2)) in bonds(poly))

_sitecolor(poly::Polyform, particle, site) =
    Roly.color(bindingrules(poly), Roly.SpeciesSiteLoc(poly.particles[particle].speciesindex, site))

"""
    twistcorrection(poly::Polyform, bond)

Return the rotation that carries the partner site of `bond` into the twist the potential is written
for.

A site with a rotational symmetry about its normal is only determined up to that symmetry, so a bond
may sit in any of `Roly.twistfreedom` twists, all of them physically the same contact. The potentials
here are written for the reference twist alone, and evaluating one at a bond that sits in another
branch puts the bonded configuration at a saddle rather than a minimum.

Throws if the sites are in no admissible twist at all, which means the contact is not a contact.
"""
function twistcorrection(poly::Polyform, bond)
    (p1, s1), (p2, s2) = bond
    D = dimension(poly)
    F = Roly.numtype(poly)
    sys = bindingrules(poly)

    b1 = bindingsite(poly.particles[p1], sys, s1)
    b2 = bindingsite(poly.particles[p2], sys, s2)
    ntwists = Roly.twistfreedom(b1, b2)
    t = Roly.twist(b1, b2)
    isnothing(t) && throw(ArgumentError(
        "sites ($p1, $s1) and ($p2, $s2) are bonded but their frames match no admissible twist " *
        "of the $ntwists this bond allows, so the bonded configuration is not a minimum"))

    return Roly.standard_rotation(F, Val(D), t, ntwists) * inv(Roly.standard_rotation(F, Val(D)))
end

"""
    map_potential(bond_potential, poly::Polyform; embed3d=false)

Build the function that gives the total bond energy of `poly` as a function of a `dtot x nparticles`
matrix of displacements from the bonded configuration.
"""
function map_potential(bond_potential, poly::Polyform; embed3d=false)
    sys = bindingrules(poly)
    bs = collect(bonds(poly))
    twists = [twistcorrection(poly, b) for b in bs]

    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    function energy_fn(ξs::AbstractMatrix)
        E = 0.0
        rottype = d == 2 ? Angle2d{eltype(ξs)} : RotXYZ{eltype(ξs)}
        for (bi, ((p1, s1), (p2, s2))) in enumerate(bs)
            particle1 = poly.particles[p1]
            particle2 = poly.particles[p2]
            spcs1 = particle1.speciesindex
            spcs2 = particle2.speciesindex

            particle_pose1 = embed3d ? Pose{3}(particle1.pose) : particle1.pose
            particle_pose2 = embed3d ? Pose{3}(particle2.pose) : particle2.pose
            site_pose1 = let sp=bindingsite(species(sys, spcs1), s1).pose
                embed3d ? Pose{3}(sp) : sp
            end
            site_pose2 = let sp=bindingsite(species(sys, spcs2), s2).pose * twists[bi]
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

"""
    hessian(bond_potential, poly; embed3d=false, perturb=false)

Hessian of [`map_potential`](@ref) at the bonded configuration, of shape
`(dtot*nparticles, dtot*nparticles)`.

`perturb=true` evaluates it a distance `√eps` away from the bonded configuration instead, for potentials
whose automatic derivatives are not finite there.
"""
function hessian(bond_potential, poly; embed3d=false, perturb=false)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    dtot = d * (d + 1) ÷ 2

    energy_fn = map_potential(bond_potential, poly; embed3d)
    x0 = zeros(Roly.numtype(poly), dtot, n)
    return _bondedhessian(energy_fn, x0; perturb)
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

"""
    EntropySolver

Supertype of the methods for turning a structure and a bond potential into a partition function.
See [`TreeLike`](@ref), [`MeanField`](@ref), [`TetheredLaplace`](@ref) and [`COMLaplace`](@ref).
"""
abstract type EntropySolver end

"""
    TetheredLaplace()

Entropy solver that expands the bond potential to second order about the bonded configuration,
holding one particle fixed to remove the global translations and rotations.
"""
struct TetheredLaplace <: EntropySolver end

"""
    COMLaplace()

Entropy solver that expands the bond potential to second order about the bonded configuration,
working in the frame of the structure's centre of mass.
"""
struct COMLaplace <: EntropySolver end

"""
    MeanField()

Entropy solver that applies mean field theory and averages over all particle tetherings.
"""
struct MeanField <: EntropySolver end

"""
    MeanFieldTopological(omega=1)
    MeanFieldCompositional(omega=1)

Entropy solver that is [`MeanField`](@ref) with every single-particle Hessian replaced by `zᵢ` copies
of the single-bond one, reducing the answer to the tree-like result times a penalty fixed by the
degree sequence alone.
"""
struct MeanFieldTopological{uniform} <: EntropySolver
    omega::Float64
end
MeanFieldTopological{uniform}() where {uniform} = MeanFieldTopological{uniform}(1.0)
MeanFieldTopological(omega::Real=1.0) = MeanFieldTopological{false}(omega)

"""
    MeanFieldCompositional(omega=1)

Entropy solver that is [`MeanField`](@ref) with every single-particle Hessian replaced by `zᵢ` copies
of the single-bond one, reducing the answer to the tree-like result times a penalty fixed by the
degree sequence alone.

`MeanFieldCompositional` substitutes the structure mean coordination number `z̄ = 2 n_b / n` for every `zᵢ`,
resulting in

    Ω_s = |SO(d)| ω^(n-1) (n / 2 n_b)^(dtot (n-1) / 2) / σ_s .

As for [`MeanFieldTopological`](@ref), `ω` is the bond volume of the bond potential, or `omega` when
there is none.

This underestimates `Ω_s`, by an amount growing with the variance of the degrees.
"""
const MeanFieldCompositional = MeanFieldTopological{true}

"""
    TreeLike(omega=1)

Entropy solver that treats every structure as a tree. A structure that really is a tree has exactly `n - 1` bonds
and the answer is exact. For non-trees, the bond volumes are averaged over the bonds. 
Given a bond potential, the bondvolumes are computed from [`bondvolume`](@ref), with no potential each 
bond contributes `omega`.
"""
struct TreeLike <: EntropySolver
    omega::Float64
end
TreeLike() = TreeLike(1.0)

"""
    EntropyModel(potential, solver=COMLaplace(); embed3d=false, perturb=false)
    EntropyModel(solver=TreeLike(); embed3d=false)

An EntropyModel is a bond potential together with a solver that computes bond volumes.

Without a potential, the solver has to be [`TreeLike`](@ref), [`MeanFieldTopological`](@ref) or
[`MeanFieldCompositional`](@ref), which then give every bond the volume `omega`.

`embed3d=true` embeds a two-dimensional structure in 3d.
`perturb=true` evaluates the Hessians of the potential a distance `√eps` away from the bonded
configuration, for potentials whose automatic derivatives are not finite there, such as ones written
with `norm(x)` at `x = 0`.
"""
struct EntropyModel{P,S<:EntropySolver}
    potential::P
    solver::S
    embed3d::Bool
    perturb::Bool
end
EntropyModel(potential, solver::EntropySolver=COMLaplace(); embed3d=false, perturb=false) =
    EntropyModel(potential, solver, embed3d, perturb)

function EntropyModel(solver::EntropySolver=TreeLike(); embed3d=false)
    solver isa Union{TreeLike,MeanFieldTopological} || throw(ArgumentError(
        "$(nameof(typeof(solver))) needs a bond potential; use `EntropyModel(potential, solver)`."))
    return EntropyModel(nothing, solver, embed3d, false)
end

"""
    entropy(model::EntropyModel, poly::Polyform)

Return the entropy `log Ω_s` of the structure `poly` under `model`, in units of `k_B`.
"""
entropy(m::EntropyModel, poly::Polyform) = _entropy(m.potential, poly, m.solver; m.embed3d, m.perturb)

# The dimension in which `m` scores the structures of `x`, a `BindingRules` or a `Polyform`: 3 whenever
# `m` embeds, and the dimension of `x` itself otherwise.
_embeddim(m::EntropyModel, x) = m.embed3d ? 3 : dimension(x)

# log(|SO(d)| / σ_s), the part of log Ω_s shared by every solver.
_logrigidvolume(poly::Polyform, d) = log(orientational_volume(d)) - log(symmetrynumber(poly))

function _entropy(::Nothing, poly::Polyform, solver::TreeLike; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    return (nparticles(poly) - 1) * log(solver.omega) + _logrigidvolume(poly, d)
end

function _entropy(bond_potential, poly::Polyform, ::TreeLike; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    logΩ = _logrigidvolume(poly, d)
    nparticles(poly) == 1 && return logΩ

    # A tree of `n` particles has exactly `n - 1` bonds, so the geometric mean is just their
    # product; anything with cycles is sees `n - 1` "typical" bonds instead.
    logω = mean(logbondvolume(bond_potential, c1, c2; d, perturb) for (c1, c2) in bondcolors(poly))
    return logΩ + (nparticles(poly) - 1) * logω
end

# The coordination number of every particle of `poly`, in one pass over its bonds.
function _degrees(poly::Polyform)
    degrees = zeros(Int, nparticles(poly))
    for bond in bonds(poly), loc in bond
        degrees[loc.particle] += 1
    end
    return degrees
end

function _entropy(bond_potential, poly::Polyform, solver::MeanFieldTopological{uniform};
                  embed3d, perturb) where {uniform}
    logtreelike = _entropy(bond_potential, poly, TreeLike(solver.omega); embed3d, perturb)
    n = nparticles(poly)
    n == 1 && return logtreelike

    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d + 1) ÷ 2
    uniform && return logtreelike - dtot * (n - 1) / 2 * log(2 * nbonds(poly) / n)
    return logtreelike - dtot * (n - 1) / (2n) * sum(log, _degrees(poly))
end

function _entropy(bond_potential, poly::Polyform, ::MeanField; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    logΩ = _logrigidvolume(poly, d)
    n == 1 && return logΩ

    _checkbondsrelaxed(bond_potential, poly)
    dtot = d * (d + 1) ÷ 2
    H = hessian(bond_potential, poly; embed3d, perturb)
    logω = 0.0
    for i in 1:n
        block = ((dtot * (i - 1)) + 1):(dtot * i)
        logω += dtot / 2 * log(2π) - logdet(Symmetric(H[block, block])) / 2
    end
    return logΩ + (n - 1) / n * logω
end

function _entropy(rsp::RigidSpringPotential, poly::Polyform, ::MeanField; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    n = nparticles(poly)
    logΩ = _logrigidvolume(poly, d)
    n == 1 && return logΩ

    _checkbondsunstressed(rsp, poly)
    logω = sum(_logomega_compatible(k, S) for (k, S) in _patchcovariances(rsp, poly, embed3d))
    return logΩ + (n - 1) / n * logω
end

function _entropy(bond_potential, poly::Polyform, ::TetheredLaplace; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d, perturb)
    λs = eigvals(H[dtot+1:end, dtot+1:end])
    S_vib = -0.5 * sum(log, λs / (2π); init=0)
    return _logrigidvolume(poly, d) + S_vib
end

function _entropy(bond_potential, poly::Polyform, ::COMLaplace; embed3d, perturb)
    d = embed3d ? 3 : dimension(poly)
    N = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    _checkbondsrelaxed(bond_potential, poly)
    H = hessian(bond_potential, poly; embed3d, perturb)
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
    return _logrigidvolume(poly, d) + log(detG) / 2 + S_vib
end
