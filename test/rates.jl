@testset "rates" begin
    # A single sphere must reproduce Stokes drags, 6πηr and 8πηr³
    r = 0.7
    θ = Crafts.frictiontensor([[0.0, 0, 0]], [r])
    @test θ ≈ θ'
    @test θ[1, 1] ≈ θ[2, 2] ≈ θ[3, 3] ≈ 6π * r
    @test θ[4, 4] ≈ θ[5, 5] ≈ θ[6, 6] ≈ 8π * r^3
    @test all(iszero, θ[1:3, 4:6]) # no translation-rotation coupling for a sphere
    Dl, Dr = Crafts.diffusionconstants([[0.0, 0, 0]], [r])
    @test Dl ≈ inv(6π * r)
    @test Dr ≈ inv(8π * r^3)

    # The two Rotne-Prager-Yamakawa branches have to meet where they change over
    for (ri, rj) in ((0.5, 0.5), (1.0, 1.0)), dir in ([1.0, 0, 0], [1, 1, 1] / sqrt(3), [0, 0, 1.0])
        R = ri + rj
        @test Crafts.rotneprager((R + 1e-9) * dir; ri, rj) ≈ Crafts.rotneprager((R - 1e-9) * dir; ri, rj) atol = 1e-8
    end
    @test Crafts.rotneprager([0.0, 0, 0]; ri=0.5, rj=0.5) ≈ Matrix(I(3) * Crafts.mobility(0.5))
    @test_throws ArgumentError Crafts.rotneprager([0.4, 0, 0]; ri=0.3, rj=0.9)

    # Friction is a tensor, so rotating the whole cluster conjugates it rather than changing it.
    Q = [0 -1 0; 1 0 0; 0 0 1.0]
    cluster, rcl = [[0.0, 0, 0], [1.1, 0, 0], [0.4, 0.9, 0]], [0.5, 0.5, 0.5]
    θc = Crafts.frictiontensor(cluster, rcl)
    θq = Crafts.frictiontensor([Q * x for x in cluster], rcl)
    @test θq ≈ [Q zeros(3, 3); zeros(3, 3) Q] * θc * [Q zeros(3, 3); zeros(3, 3) Q]'
    @test θc ≈ θc'

    # Friction cannot depend on how the particles are numbered.
    cl, rcl = [[0.0, 0, 0], [1.1, 0, 0], [0.4, 0.9, 0]], [0.5, 0.5, 0.5]
    @test Crafts.frictiontensor(cl, rcl) ≈ Crafts.frictiontensor(cl[[3, 1, 2]], rcl[[3, 1, 2]])
    @test eltype(Crafts.frictiontensor(Vector{Float32}.(cl), Float32.(rcl))) == Float32

    # The center of diffusion is by definition the point where the coupling block is symmetric.
    for (xs, rs) in ((cl, rcl), ([[-1.0, 0, 0], [1.0, 0, 0]], [0.3, 0.9]))
        D = Crafts.diffusiontensor(xs, rs)
        @test D[1:3, 4:6] ≈ D[1:3, 4:6]' atol = 1e-12 * maximum(abs, D)
    end

    # Widening the caps to cover the spheres recovers the isotropic Smoluchowski limit.
    @test Crafts.orientedbindingrate(; D=0.4, Dr1=0.1, δ1=π, Dr2=0.1, δ2=π, R=2.0) ≈
          Crafts.smoluchowskirate(; D=0.4, R=2.0)
    # and nothing oriented can ever beat it
    @test all(Crafts.orientedbindingrate(; D=0.4, Dr1=0.1, δ1=d, Dr2=0.1, δ2=d, R=2.0) <=
              Crafts.smoluchowskirate(; D=0.4, R=2.0) for d in 0.05:0.05:π)

    # A pair mobility may never exceed the self mobility
    for R in 0.05:0.05:2.0
        @test opnorm(Crafts.rotneprager([R, 0, 0]; ri=0.5, rj=0.5)) <= Crafts.mobility(0.5) + 1e-12
    end

    # Regression: spheres this close are indefinite under the far-field form alone, and used to give
    # negative diffusion constants and negative binding rates.
    overlapping = [[0.0, 0, 0], [0.5774, 0, 0], [0.2887, 0.5, 0], [0.8661, 0.5, 0]]
    rov = fill(0.5774, 4)
    @test minimum(norm(overlapping[a] - overlapping[b]) for a in 1:4 for b in (a+1):4) < 2rov[1]
    @test isposdef(Symmetric(Crafts.frictiontensor(overlapping, rov)))
    @test all(>(0), Crafts.diffusionconstants(overlapping, rov))

    # The center of diffusion sits at the midpoint by symmetry, and shifts to the larger sphere.
    dimer = [[-1.0, 0, 0], [1.0, 0, 0]]
    @test Crafts.centerofdiffusion(dimer, [0.5, 0.5]) ≈ [0, 0, 0] atol = 1e-12
    @test Crafts.centerofdiffusion(dimer, [0.3, 0.9])[1] > 0

    # Adding particles can only slow a cluster down.
    chain(n) = Crafts.diffusionconstants([[1.0i, 0, 0] for i in 1:n], fill(0.5, n))
    Ds = [chain(n) for n in 1:4]
    @test issorted(first.(Ds); rev=true)
    @test issorted(last.(Ds); rev=true)

    # The hydrodynamic radius is the inradius of each polygon, so bonded particles meet at contact
    # rather than interpenetrating as they would at the circumradius.
    @test Crafts.sitedistance(UnitTriangle) ≈ 1 / (2sqrt(3))
    @test Crafts.sitedistance(UnitSquare) ≈ 1 / 2
    @test Crafts.sitedistance(UnitHexagon) ≈ sqrt(3) / 2
    @test all(s -> Crafts.sitedistance(s) < bounding_radius(s), (UnitTriangle, UnitSquare, UnitHexagon))

    systems = (BindingRules([1 3 2 3; 2 1 4 1; 2 2 3 2; 3 1 4 1], UnitTriangle),
               BindingRules([1 1 2 3; 2 1 3 3; 1 2 4 4; 2 2 5 4], UnitSquare),
               BindingRules([1 1 2 4; 2 2 3 5; 3 3 1 6], UnitHexagon))
    for rules in systems
        asys = AssemblySystem(rules, EntropyModel(TreeApproximation()); maxsize=5, verbose=false)
        net = ReactionNetwork(asys)

        for rxn in net
            p = product(rxn)
            xs, rs = Crafts.particlepositions(p), Crafts.particleradii(p)

            # every fragment must be a physically realisable rigid body
            for h in halves(rxn)
                idx = Crafts.particlesat(p, h)
                @test !isempty(idx)
                @test isposdef(Symmetric(Crafts.frictiontensor(xs[idx], rs[idx])))
            end
            # the halves partition the particles
            @test sort(vcat(Crafts.particlesat(p, halves(rxn)[1]),
                            Crafts.particlesat(p, halves(rxn)[2]))) == 1:nparticles(p)

            # the broken bond's two sites sit on top of one another
            for e in cut(rxn)
                @test Crafts.siteposition(p, e.src) ≈ Crafts.siteposition(p, e.dst) atol = 1e-8
            end
        end

        # `denom` has a inv(1 - Λ) term. test that Λ cannot reach 1 (it is actually bounded by 1/√2)
        Λ(δ, ξ) = sin(δ / 2)^2 * (ξ + cot(δ / 2)) / (ξ + sin(δ / 2) * cos(δ / 2))
        @test maximum(Λ(δ, ξ) for δ in range(1e-6, π / 2; length=500), ξ in range(1 / sqrt(2), 50; length=500)) ≤ 1 / sqrt(2)

        rates = [orientedbindingrate(rxn) for rxn in net]
        @test all(x -> isfinite(x) && x > 0, rates)
        @test all(x -> isfinite(x) && x > 0, [smoluchowskirate(rxn) for rxn in net])

        # a wider reactive patch can only help
        @test all(orientedbindingrate(rxn; siteradius=1) > orientedbindingrate(rxn; siteradius=0.5) for rxn in net)

        # usable directly as a kernel, and doing so leaves the network detailed balanced
        rate!(net; fwdkernel=orientedbindingrate)
        @test fwdrates(net) ≈ rates
        @test isdetailedbalanced(net)
        @test all(net.active)
    end
end;
