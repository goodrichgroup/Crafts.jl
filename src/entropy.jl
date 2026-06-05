function map_potential(bond_potential, poly::Polyform)
    sys = assemblysystem(poly)
    bs = collect(bonds(poly))

    d = dimension(poly)
    dtot = d * (d+1) ÷ 2

    function energy_fn(ξs::AbstractMatrix)
        E = 0.0
        rottype = d == 2 ? Angle2d{eltype(ξs)} : RotXYZ{eltype(ξs)}
        for ((p1, s1), (p2, s2)) in bs
            particle1 = poly.particles[p1]
            particle2 = poly.particles[p2]
            spcs1 = particle1.species_index
            spcs2 = particle2.species_index

            # transform coordinates such that minimum is at zero.
            pose1 = Pose(ξs[1:d, p1] + particle1.pose.x, rottype(ξs[d+1:dtot, p1]...) * particle1.pose.psi)
            pose2 = Pose(ξs[1:d, p2] + particle2.pose.x, rottype(ξs[d+1:dtot, p2]...) * particle2.pose.psi)

            x1 = pose1 * bindingsites(species(sys, spcs1), s1).pose
            x2 = pose2 * bindingsites(species(sys, spcs2), s2).pose

            E += bond_potential(x1, x2)
        end
        return E
    end

    return energy_fn
end

function hessian(bond_potential, poly)
    d = dimension(poly)
    n = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    energy_fn = map_potential(bond_potential, poly)
    x0 = zeros(numtype(poly), dtot, n)

    return ForwardDiff.hessian(energy_fn, x0 .+ sqrt(eps(eltype(x0))))
end

function orientational_volume(d::Integer)
    return d == 2 ? 2π : d == 3 ? 8π^2 : throw(ArgumentError("orientational volume requires d=2 or d=3"))
end

function inertia_tensor(poly::Polyform)
    d = dimension(poly)
    com = sum(let x=p.pose.x; x end for p in poly.particles) / nparticles(poly)
    return sum(let x=p.pose.x - com; norm(x)^2 * I(d) - x * x' end for p in poly.particles)
end

function rotational_inertia(poly::Polyform{2})
    return tr(inertia_tensor(poly)) * I(1)
end

function rotational_inertia(poly::Polyform{3})
   return inertia_tensor(poly)
end

abstract type EntropyApproximator end
struct TetheredLaplace <: EntropyApproximator end
struct COMLaplace <: EntropyApproximator end


function entropy(bond_potential, poly::Polyform; method::EntropyApproximator=TetheredLaplace())
    return _entropy(bond_potential, poly, method)
end

function _entropy(bond_potential, poly::Polyform, ::TetheredLaplace)
    d = dimension(poly)
    dtot = d * (d+1) ÷ 2

    H = hessian(bond_potential, poly)
    λs = eigvals(H[dtot+1:end, dtot+1:end])
    S_vib = -0.5 * sum(log, λs / (2π); init=0)
    return orientational_volume(d) * exp(S_vib) / symmetrynumber(poly)
end

function _entropy(bond_potential, poly::Polyform, ::COMLaplace)
    d = dimension(poly)
    N = nparticles(poly)
    dtot = d * (d+1) ÷ 2

    H = hessian(bond_potential, poly)
    λs = eigvals(H)
    S_vib = -0.5 * sum(log, λ / (2π) for λ in @view λs[dtot+1:end]; init=0)

    inertia = inertia_tensor(poly)
    if d == 2
        detG = N^d * (tr(inertia) + N)
    elseif d == 3
        detG = N^d * det(inertia + N * I(dtot - d))
    else
        throw(ArgumentError("dimension must be 2 or 3"))
    end
    return orientational_volume(d) * sqrt(detG) * exp(S_vib) / symmetrynumber(poly)
end