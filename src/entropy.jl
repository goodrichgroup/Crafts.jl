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

            E += bond_potential(x1, x2)
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

function entropy(bond_potential, poly::Polyform; method::EntropySolver=TetheredLaplace(), embed3d=false)
    return _entropy(bond_potential, poly, method; embed3d)
end

function _entropy(bond_potential, poly::Polyform, ::TetheredLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    dtot = d * (d+1) ÷ 2

    H = hessian(bond_potential, poly; embed3d)
    λs = eigvals(H[dtot+1:end, dtot+1:end])
    S_vib = -0.5 * sum(log, λs / (2π); init=0)
    return orientational_volume(d) * exp(S_vib) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::COMLaplace; embed3d)
    d = embed3d ? 3 : dimension(poly)
    N = nparticles(poly)
    dtot = d * (d+1) ÷ 2

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
