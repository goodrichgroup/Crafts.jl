mobility(r) = inv(6π * r)

"""
    rotneprager(Δx; ri=1, rj=1)

Rotne-Prager-Yamakawa pair mobility of two spheres separated by `Δx`.

Overlapping spheres take Yamakawa's near-field branch. The far-field form is not positive definite
there: it returns pair mobilities exceeding the self mobility, which makes the grand mobility matrix
indefinite and diffusion constants come out negative. The two branches agree at contact.
"""
function rotneprager(Δx; ri=1, rj=1)
    R = norm(Δx)
    a = (ri + rj) / 2
    iszero(R) && return Matrix{eltype(Δx)}(I(3) / (6π * a))

    P = Δx * Δx' / R^2
    R >= ri + rj && return (I + P + (ri^2 + rj^2) / R^2 * (I / 3 - P)) / (8π * R)

    ri ≈ rj || throw(ArgumentError("cannot construct RPY tensor for overlapping spheres of unequal size"))
    return ((1 - 9R / (32a)) * I + (3R / (32a)) * P) / (6π * a)
end

"""
    frictiontensor(xs, rs)

The `6 x 6` grand friction tensor of spheres of radii `rs` at positions `xs`, held rigidly together.
Blocks are translation-translation, translation-rotation and rotation-rotation, about the origin.
"""
function frictiontensor(xs, rs)
    n = length(xs)
    V = 4 / 3 * π * sum(x -> x^3, rs)

    T = float(promote_type(eltype(eltype(xs)), eltype(rs)))
    B = zeros(T, 3n, 3n)
    K = zeros(T, 3n, 6)
    for i in 1:n
        B[1+3(i-1):3i, 1+3(i-1):3i] .= I(3) * mobility(rs[i])
        for j in (i+1):n
            t = rotneprager(xs[j] - xs[i]; ri=rs[i], rj=rs[j])
            B[1+3(i-1):3i, 1+3(j-1):3j] .= t
            B[1+3(j-1):3j, 1+3(i-1):3i] .= t
        end

        x, y, z = xs[i]
        K[1+3(i-1):3i, 1:3] .= I(3)
        K[1+3(i-1):3i, 4:6] .= [0 z -y; -z 0 x; y -x 0]
    end

    #   C = inv(B),   E_i = I,   V_i = U(i)',   K = [E V]
    #
    #   θt        =  Σ_ij C_ij              =  E' C E
    #   θtr       =  Σ_ij U(i) C_ij         =  V' C E
    #   θrr_uncor = -Σ_ij U(i) C_ij U(j)    =  V' C V
    #
    #   [θt θtr'; θtr θrr_uncor]            =  K' C K  =  K' * (B \ K)

    # `cholesky` refuses an indefinite mobility, extra protection against out-of-range Rotne-Prager
    θ = K' * (cholesky(Symmetric(B)) \ K)
    θ[4:6, 4:6] += 6V * I
    return θ
end

"""
    centerofdiffusion(xs, rs)

The point about which the translation-rotation coupling of the cluster is symmetric.
"""
function centerofdiffusion(xs, rs)
    D0 = inv(frictiontensor(xs, rs))

    dtr = @view D0[4:6, 1:3]
    drr = @view D0[4:6, 4:6]

    A = tr(drr) * I - drr
    b = [dtr[2, 3] - dtr[3, 2], dtr[3, 1] - dtr[1, 3], dtr[1, 2] - dtr[2, 1]]

    r_od = A \ b
    return r_od
end

"""
    diffusiontensor(xs, rs; cod=nothing)

The `6 x 6` diffusion tensor of the cluster about its center of diffusion.
"""
function diffusiontensor(xs, rs; cod=nothing)
    if isnothing(cod)
        cod = centerofdiffusion(xs, rs)
    end
    D = inv(frictiontensor(xs .- Ref(cod), rs))
    return D
end

"""
    diffusionconstants(xs, rs; cod=nothing)
    diffusionconstants(p::Polyform; cod=nothing)

The orientationally averaged translational and rotational diffusion constants, `(Dlin, Drot)`.
"""
diffusionconstants(xs, rs; kwargs...) = let d = diffusiontensor(xs, rs; kwargs...)
    (d[1, 1] + d[2, 2] + d[3, 3]) / 3, (d[4, 4] + d[5, 5] + d[6, 6]) / 3
end

# embed 2d coordinates in 3d
_pad3(x::AbstractVector) = SVector{3}(x[1], x[2], length(x) > 2 ? x[3] : zero(eltype(x)))

particlepositions(p::Polyform) = [_pad3(part.pose.x) for part in p.particles]

# we use this for the hydrodynamic radius of a species
sitedistance(ps::ParticleSpecies) = mean(norm(bindingsites(ps, i).pose.x) for i in 1:nsites(ps))

function particleradii(p::Polyform)
    sys = bindingrules(p)
    return [sitedistance(species(sys, Roly.species_index(part))) for part in p.particles]
end

siteposition(p::Polyform, v::Integer) = let (i, j) = Roly._vertex_to_particle_site(p, v)
    _pad3(bindingsites(p.particles[i], bindingrules(p), j).pose.x)
end

# `halves` indexes binding sites, several of which belong to the same particle.
particlesat(p::Polyform, vs) = unique!([Roly._vertex_to_particle_site(p, v)[1] for v in vs])

function diffusionconstants(p::Polyform; kwargs...)
    return diffusionconstants(particlepositions(p), particleradii(p); kwargs...)
end

##### Diffusion-limited association rates #####

"""
    capfraction(δ)

Fraction of a sphere covered by a cap of half-angle `δ`.
"""
capfraction(δ) = sin(δ / 2)^2

"""
    orientedbindingrate(; D, Dr1, δ1, Dr2, δ2, R)
    orientedbindingrate(rxn::Reaction; siteradius=0.5)

Diffusion-limited association rate when binding also requires the two partners to be correctly
oriented: reactive caps of half-angle `δ1`, `δ2` on spheres of separation `R`, with relative
translational diffusion `D` and rotational diffusion `Dr1`, `Dr2`.

Reduces to [`smoluchowskirate`](@ref) as the caps grow to cover their spheres.
"""
function orientedbindingrate(; D, Dr1, δ1, Dr2, δ2, R)
    # for strongly non-convex bodies, it may be the case that R1 + R2 != R
    F1 = capfraction(δ1)
    F2 = capfraction(δ2)

    ξ1 = sqrt((1 + Dr1 * R^2 / D) / 2)
    ξ2 = sqrt((1 + Dr2 * R^2 / D) / 2)

    Λ1 = F1 * (ξ1 + cot(δ1 / 2)) / (ξ1 + sin(δ1 / 2) * cos(δ1 / 2))
    Λ2 = F2 * (ξ2 + cot(δ2 / 2)) / (ξ2 + sin(δ2 / 2) * cos(δ2 / 2))

    denom = Λ1 * Λ2 +
            inv(inv(1 - Λ1) * inv(1 - Λ2) + inv(1 - Λ1) * inv(Λ2 - F2) + inv(1 - Λ2) * inv(Λ1 - F1))
    k = 4π * R * D * F1 * F2 / denom
    return k
end

"""
    orientedbindingrate(rxn::Reaction; siteradius=0.5)

Diffusion-limited association rate when binding also requires the two partners to be correctly
oriented. Takes relevant parameters from the geometry reacting particles. `siteradius` determines the size of the
binding patches.

Usable directly as a kernel for [`rate!`](@ref).
"""
function orientedbindingrate(rxn::Reaction; siteradius=0.5)
    p = product(rxn)
    xs, rs = particlepositions(p), particleradii(p)

    # Center of binding: the sites being broken sit on top of one another, so either endpoint would work
    cob = mean(siteposition(p, e.src) for e in cut(rxn))

    vs1, vs2 = halves(rxn)
    idxs1, idxs2 = particlesat(p, vs1), particlesat(p, vs2)

    xs1, rs1 = xs[idxs1], rs[idxs1]
    xs2, rs2 = xs[idxs2], rs[idxs2]

    cod1 = centerofdiffusion(xs1, rs1)
    cod2 = centerofdiffusion(xs2, rs2)

    Dl1, Dr1 = diffusionconstants(xs1, rs1; cod=cod1)
    Dl2, Dr2 = diffusionconstants(xs2, rs2; cod=cod2)

    R = norm(cod2 - cod1)
    r1 = norm(cob - cod1)
    r2 = norm(cob - cod2)

    R1 = abs((r1^2 - r2^2 + R^2) / (2R))
    R2 = abs((r2^2 - r1^2 + R^2) / (2R))

    δ1 = atan(siteradius / R1)
    δ2 = atan(siteradius / R2)
    D = Dl1 + Dl2
    return orientedbindingrate(; D, Dr1, δ1, Dr2, δ2, R)
end

"""
    smoluchowskirate(; D, R)

Diffusion-limited association rate of two spheres that react on contact whatever their orientation,
`4π D R`.
"""
smoluchowskirate(; D, R) = 4π * D * R

"""
    smoluchowskirate(rxn::Reaction)

Diffusion-limited association rate of two spheres that react on contact whatever their orientation,
`4π D R`. Usable directly as a kernel for [`rate!`](@ref).
"""
function smoluchowskirate(rxn::Reaction)
    p = product(rxn)
    xs, rs = particlepositions(p), particleradii(p)

    vs1, vs2 = halves(rxn)
    idxs1, idxs2 = particlesat(p, vs1), particlesat(p, vs2)

    xs1, rs1 = xs[idxs1], rs[idxs1]
    xs2, rs2 = xs[idxs2], rs[idxs2]

    cod1 = centerofdiffusion(xs1, rs1)
    cod2 = centerofdiffusion(xs2, rs2)

    Dl1, _ = diffusionconstants(xs1, rs1; cod=cod1)
    Dl2, _ = diffusionconstants(xs2, rs2; cod=cod2)

    R = norm(cod2 - cod1)
    D = Dl1 + Dl2
    return smoluchowskirate(; D, R)
end
