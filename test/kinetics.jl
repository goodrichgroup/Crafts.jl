@testset "kinetics" begin
    dimer = BindingRules([1 1 2 3], UnitTriangle)   # 2 species, 1 bond: A, B, AB
    rules = BindingRules(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangle)   # 16 structures, largest has 5 particles
    model = EntropyModel(TreeLike())

    geomean(xs) = exp(sum(log, xs) / length(xs))

    asys = AssemblySystem(dimer, model; verbose=false)
    net = ReactionNetwork(asys)
    ξ = [-2.0, -2.5, 3.0]

    # The single reaction A + B <-> AB, so the scale is just its detailed-balance factor
    _, kscale = Crafts._make_odefunction_and_ratescale(net, ξ)
    ρs = densities(asys, ξ)
    @test kscale ≈ ρs[1] * ρs[2] / ρs[3]

    # chemical potentials cancel out of the detailed-balance factor
    for shift in (5.0, -3.0)
        _, shifted = Crafts._make_odefunction_and_ratescale(net, [ξ[1:2] .+ shift; ξ[3]])
        @test shifted ≈ kscale
    end
    _, colder = Crafts._make_odefunction_and_ratescale(net, [ξ[1:2]; 5.0])
    @test !(colder ≈ kscale)

    asys16 = AssemblySystem(rules, model; verbose=false)
    ξ16 = [fill(-3.0, 4); fill(2.0, 4)]

    # A reaction can be active on its forward rate alone; such rates drop out of the mean
    irreversible(rxn) = nparticles(product(rxn)) == 2 ? 0.0 : 1.0
    net16 = ReactionNetwork(asys16; bwdkernel=irreversible)
    @test all(net16.active)

    logρs = logdensities(asys16, ξ16)
    factors = [Crafts._copy_active_bwdrates(net16)[r] * exp(logρs[i] + logρs[j] - logρs[k])
               for (r, (i, j, k)) in enumerate(Crafts._copy_active_reactions(net16))]
    @test any(iszero, factors)
    _, kscale16 = Crafts._make_odefunction_and_ratescale(net16, ξ16)
    @test kscale16 ≈ geomean(filter(>(0), factors))

    # with no dissociation at all the forward rates set the scale instead
    onesided = ReactionNetwork(asys16; fwdkernel=rxn -> inv(length(cut(rxn))), bwdkernel=Returns(0.0))
    @test all(onesided.active)
    stepo, kscaleo = Crafts._make_odefunction_and_ratescale(onesided, ξ16)
    @test kscaleo ≈ geomean(Crafts._copy_active_fwdrates(onesided))

    @test length(reactions(net16)) == nreactions(net16)
    @test reactions(net16)[3] == indices(net16[3])
    @test reactions(net16) == map(indices, collect(net16))
    @test length(Crafts._copy_active_reactions(net16)) == count(net16.active)
    reactions(net16), fwdrates(net16), bwdrates(net16)
    @test (@allocated reactions(net16)) == 0
    @test (@allocated fwdrates(net16)) == 0
    @test (@allocated bwdrates(net16)) == 0

    # irreversible aggregation only ever consumes structures, but still conserves particles
    _, uso = simulate_kinetics(onesided, ξ16; T=1e3)
    @test all(≥(-1e-8), uso)
    @test all(≤(1e-8), diff(vec(sum(uso; dims=1))))
    ϕo = @view(compositionmatrix(asys16)[:, 1:4])' * uso
    @test all(≈(ϕo[:, 1]), eachcol(ϕo))

    # test invariance of steady state to kernels
    kernels = (Returns(1.0),
               rxn -> 1.0 + nparticles(product(rxn)),
               rxn -> 2.0^sum(bondcounts(rxn)),
               rxn -> inv(length(cut(rxn))))
    for a in (asys, asys16)
        ξa = [fill(-3.0, nspecies(a)); fill(2.0, nbonds(a))]
        ρeq = densities(a, ξa)
        du = similar(ρeq)

        for fwdkernel in kernels
            n = ReactionNetwork(a; fwdkernel)
            step, scale = Crafts._make_odefunction_and_ratescale(n, ξa)

            @test isdetailedbalanced(n)
            step(du, ρeq / scale, nothing, 0.0)
            @test maximum(abs, du) < 1e-12
        end

        # unequal kernels break equilibrium, which is exactly what `isdetailedbalanced` reports
        driven = ReactionNetwork(a; fwdkernel=Returns(1.0), bwdkernel=Returns(4.0))
        @test !isdetailedbalanced(driven)
        step, scale = Crafts._make_odefunction_and_ratescale(driven, ξa)
        step(du, ρeq / scale, nothing, 0.0)
        @test maximum(abs, du) > 1e-3

        # an inactive reaction cannot unbalance anything, since it contributes nothing either way
        @test isdetailedbalanced(ReactionNetwork(a; fwdkernel=Returns(0.0), bwdkernel=Returns(0.0)))
    end

    ts, us = simulate_kinetics(net, ξ; T=50.0)
    @test length(ts) == size(us, 2)
    @test size(us, 1) == nstructures(asys)
    @test ts[1] == 0 && ts[end] ≈ 50.0
    @test us[:, 1] ≈ [particledensities(asys, ξ); 0.0]
    @test us[:, end] ≈ ρs rtol = 1e-6
    @test length(first(simulate_kinetics(net, ξ; T=50.0, saveat=0:5:50))) == 11

    # particle densities are conserved along the whole trajectory, not just at its ends
    N = @view compositionmatrix(asys)[:, 1:nspecies(asys)]
    ϕts = N' * us
    @test all(≈(ϕts[:, 1]), eachcol(ϕts))

    # entering through concentrations must reproduce the concentrations asked for
    ϕs, εs = [0.4, 0.3], [3.0]
    _, usϕ = simulate_kinetics(net, ϕs, εs; T=50.0)
    @test N' * usϕ[:, end] ≈ ϕs rtol = 1e-6

    # A system of nothing but monomers and dimers has one bimolecular reaction, so eliminating the
    # conserved particle densities leaves a scalar Riccati equation with an elementary solution.
    # That gives the solver an exact trajectory to be checked against rather than only a fixed point.

    # A + B ⟷ AB:  dρ/dt = kf(a - ρ)(b - ρ) - kb ρ = kf (ρ - ρ₊)(ρ - ρ₋),  ρ = ρ_AB, ρ(0) = 0
    let n = ReactionNetwork(asys; fwdkernel=Returns(2.5)), ξd = [-2.0, -2.5, 3.0]
        ρeq = densities(asys, ξd)
        kf = fwdrates(n)[1]
        kb = bwdrates(n)[1] * ρeq[1] * ρeq[2] / ρeq[3] # the ODE's backward rate carries detailed balance
        a, b = particledensities(asys, ξd)

        # the discriminant of that quadratic, so ρ± are its roots and ρ₊ - ρ₋ = disc/kf the rate
        disc = sqrt((kf * (a + b) + kb)^2 - 4kf^2 * a * b)
        ρp, ρm = (kf * (a + b) + kb + disc) / 2kf, (kf * (a + b) + kb - disc) / 2kf
        ρAB(t) = let E = ρm / ρp * exp(-kf * (ρp - ρm) * t)
            (ρm - E * ρp) / (1 - E)
        end

        tds = [0.0, 0.05, 0.2, 0.5, 1.0, 3.0, 10.0]
        _, ud = simulate_kinetics(n, ξd; T=10.0, saveat=tds, reltol=1e-12, abstol=1e-14)
        @test ud[3, :] ≈ ρAB.(tds) rtol = 1e-10
        @test ud[1, :] ≈ a .- ρAB.(tds) rtol = 1e-10
        @test ud[2, :] ≈ b .- ρAB.(tds) rtol = 1e-10
        # the root it runs into is the equilibrium, and the rate it approaches it at is the only mode
        @test ρm ≈ ρeq[3]
        @test correlationtime(n, ξd) ≈ inv(kf * (ρp - ρm))
    end

    # 2A ⟷ A₂, where `i == j` makes the reaction land on `du[i]` twice:
    # dρ/dt = 2(kb ρ_A₂ - kf ρ²) = kb(c - ρ) - 2kf ρ² = -2kf (ρ - ρ₊)(ρ - ρ₋),  ρ = ρ_A, ρ(0) = c
    let homo = AssemblySystem(BindingRules([1 1 1 1], UnitTriangle), model; verbose=false), ξh = [-1.5, 2.0]
        @test nparticles.(structures(homo)) == [1, 2]
        n = ReactionNetwork(homo; fwdkernel=Returns(1.7))
        @test reactions(n) == [(1, 1, 2)]

        ρeq = densities(homo, ξh)
        kf = fwdrates(n)[1]
        kb = bwdrates(n)[1] * ρeq[1]^2 / ρeq[2]
        c = particledensities(homo, ξh)[1]

        # the discriminant of that quadratic, so ρ± are its roots and ρ₊ - ρ₋ = disc/2kf the rate
        disc = sqrt(kb^2 + 8kf * kb * c)
        ρp, ρm = (-kb + disc) / 4kf, (-kb - disc) / 4kf
        E0 = (c - ρp) / (c - ρm)
        ρA(t) = let E = E0 * exp(-2kf * (ρp - ρm) * t)
            (ρp - ρm * E) / (1 - E)
        end

        ths = [0.0, 0.02, 0.1, 0.3, 1.0, 5.0]
        _, uh = simulate_kinetics(n, ξh; T=5.0, saveat=ths, reltol=1e-12, abstol=1e-14)
        @test uh[1, :] ≈ ρA.(ths) rtol = 1e-10
        @test uh[2, :] ≈ (c .- ρA.(ths)) ./ 2 rtol = 1e-10
        @test ρp ≈ ρeq[1]
        @test correlationtime(n, ξh) ≈ inv(2kf * (ρp - ρm))
    end

    # The rates do not depend on ξ's implied particle densities, so initial densities that disagree with it
    # simply relax to the equilibrium implied by their own
    ξϕ = topotentials(asys, ϕs, εs)
    _, usρ0 = simulate_kinetics(net, ξϕ; T=50.0, initial_densities=[2ϕs; 0.0])
    @test usρ0[:, end] ≈ densities(asys, topotentials(asys, 2ϕs, εs)) rtol = 1e-6
    @test !isapprox(usρ0[:, end], densities(asys, ξϕ); rtol=1e-3)
end;
