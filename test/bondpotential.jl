@testset "bondpotential" begin
    rng = MersenneTwister(20260806)
    randpose(::Val{2}) = Pose(SVector{2}(randn(rng, 2)), Angle2d(2π * rand(rng)))
    randpose(::Val{3}) = Pose(SVector{3}(randn(rng, 3)), RotXYZ(2π .* rand(rng, 3)...))
    # every cloud is used through its spread about its own centroid, so compare against centred ones
    centred(D, n) = let A = randn(rng, D, n)
        A .- mean(A; dims=2)
    end

    # The stored covariances are only a rewriting of the raw spring sum, so the two must agree.
    springs(A, B, ks) = (p1, p2) -> sum(ks[i] / 2 * sum(abs2, (p1 * A[:, i]) - (p2 * B[:, i]))
                                        for i in axes(A, 2))
    agrees(f, g, D) = all(1:50) do _
        p1, p2 = randpose(Val(D)), randpose(Val(D))
        isapprox(f(p1, p2), g(p1, p2); atol=1e-12, rtol=1e-12)
    end

    for D in (2, 3), n in (3, 5)
        A, B = centred(D, n), centred(D, n)
        k, ks = 3.7, rand(rng, n) .+ 0.5

        two = RigidSpringPotential([A, B]; k)
        @test agrees((p1, p2) -> two(p1, p2, 1, 2), springs(A, B, fill(k / n, n)), D)
        # the two sites may be handed over in either order
        @test agrees((p1, p2) -> two(p1, p2, 1, 2), (p1, p2) -> two(p2, p1, 2, 1), D)
        @test two.Sij[2, 1] ≈ two.Sij[1, 2]'

        # per-patch stiffnesses also move the centroid the potential centres on
        w = ks ./ sum(ks)
        weighted = RigidSpringPotential([A, B]; k=ks)
        @test agrees((p1, p2) -> weighted(p1, p2, 1, 2), springs(A .- A * w, B .- B * w, ks), D)
        @test weighted.k[1, 2] ≈ sum(ks)

        # a single cloud is paired with the rotated copy that lets every spring relax at contact
        M = Matrix(Roly.standard_rotation(Float64, Val(D)))
        @test det(M) ≈ 1   # the partner is reachable by a rigid motion, so a rotation, not a mirror
        one_ = RigidSpringPotential(A; k)
        @test agrees(one_, springs(A, M * A, fill(k / n, n)), D)
        p = randpose(Val(D))
        @test one_(p, Pose(p.x, p.psi * inv(Roly.standard_rotation(Float64, Val(D))))) ≈ 0 atol = 1e-12
        @test strainenergy(one_) ≈ 0 atol = 1e-12

        # only where a cloud sits relative to its own centroid matters, not where the centroid is
        @test agrees(RigidSpringPotential(A .+ randn(rng, D); k), one_, D)
    end

    # Only the relative pose enters, so rigidly moving both sites cannot change the energy.
    let rsp = RigidSpringPotential(centred(2, 4); k=2.0)
        p1, p2, g = randpose(Val(2)), randpose(Val(2)), randpose(Val(2))
        @test rsp(g * p1, g * p2) ≈ rsp(p1, p2)
    end

    # A cloud that leaves the bond free to rotate has no stiff limit, so it is refused outright.
    @test_throws ArgumentError RigidSpringPotential(reshape([0.3, -0.2], 2, 1))       # 2d, one patch
    @test_throws ArgumentError RigidSpringPotential([0.0 1 2; 0 0 0; 0 0.5 1])        # 3d, collinear
    @test_throws ArgumentError RigidSpringPotential(0.0)                              # σ = 0
    @test_throws ArgumentError RigidSpringPotential(1.0, 0.0)
    @test RigidSpringPotential([0.0 1 0; 0 0 1; 0 0 0]) isa RigidSpringPotential      # 3d, planar
    @test_throws ArgumentError RigidSpringPotential([centred(2, 3), reshape([0.3, -0.2], 2, 1)])

    # Springs that cannot all relax at once leave strain behind, and that strain suppresses ω.
    let A = centred(2, 4), B = centred(2, 4), k = 40.0
        strained = RigidSpringPotential([A, B]; k)
        @test strainenergy(strained, 1, 2) > 0
        @test logbondvolume(strained, 1, 2) < logbondvolume(RigidSpringPotential(A; k))
        # strain scales with the stiffness, since it is an energy stored in the springs
        @test strainenergy(RigidSpringPotential([A, B]; k=10k), 1, 2) ≈ 10 * strainenergy(strained, 1, 2)
    end

    # For compatible patches the exponent cancels and ω = √(2π/k)^dtot / ∏_{μ<ν} √(λμ + λν),
    # with λ the eigenvalues of S. Both uniform forms have one zero mode along the site normal.
    for (σ, k) in ((1.0, 30.0), (0.4, 250.0))
        t = σ^2 / 12
        @test bondvolume(RigidSpringPotential(σ; k)) ≈ sqrt(2π / k)^3 / sqrt(t)
    end
    for (σx, σy, k) in ((1.0, 1.0, 30.0), (0.4, 1.3, 250.0))
        tx, ty = σx^2 / 12, σy^2 / 12
        @test bondvolume(RigidSpringPotential(σx, σy; k)) ≈
              sqrt(2π / k)^6 / sqrt(tx * ty * (tx + ty))
    end
    # so ω scales as k^(-dtot/2)
    @test bondvolume(RigidSpringPotential(1.0; k=1.0)) /
          bondvolume(RigidSpringPotential(1.0; k=100.0)) ≈ 100^(3 / 2)
    @test bondvolume(RigidSpringPotential(1.0, 1.0; k=1.0)) /
          bondvolume(RigidSpringPotential(1.0, 1.0; k=100.0)) ≈ 100^(6 / 2)

    # Potentials that know nothing about bond volumes have to say so rather than guess.
    @test_throws ArgumentError bondvolume((p1, p2) -> 0.0)

    # A cloud paired with an unrotated copy of itself relaxes when the two sites coincide rather
    # than when they face each other. That is a legitimate model, and `strainenergy` stays zero
    # since the springs do relax somewhere, but the bond is not at rest where Roly puts it.
    let A = centred(2, 4), k = 10.0
        elsewhere = RigidSpringPotential([A]; k)
        here = RigidSpringPotential(A; k)
        p = randpose(Val(2))
        contact = Pose(p.x, p.psi * inv(Roly.standard_rotation(Float64, Val(2))))

        @test strainenergy(elsewhere) ≈ 0 atol = 1e-12
        @test contactexcess(elsewhere) ≈ elsewhere(p, contact) > 0
        @test contactexcess(here) ≈ 0 atol = 1e-12
        @test contactexcess(elsewhere) ≈ 2k * tr(elsewhere.Scc[1])
        # a potential we cannot inspect is taken at its word
        @test contactexcess((p1, p2) -> 0.0) == 0
    end

    # A potential that cannot describe the rules it is handed is caught when the system is built,
    # rather than as an index or method error from inside a partition function.
    let rules2d = BindingRules([1 1 2 3], UnitTriangle)
        @test_throws ArgumentError AssemblySystem(rules2d, EntropyModel(RigidSpringPotential(1.0, 1.0),
                                                                       TreeLike()))
        @test_throws ArgumentError AssemblySystem(rules2d, EntropyModel(RigidSpringPotential(1.0),
                                                                       TreeLike();
                                                                       embed3d=true))
        @test_throws ArgumentError AssemblySystem(rules2d,
                                                  EntropyModel(RigidSpringPotential([centred(2, 3)]),
                                                               TreeLike()))
        # a color-independent potential covers any number of colors, and closures are not inspected
        @test AssemblySystem(rules2d, EntropyModel(RigidSpringPotential(1.0), TreeLike())) isa
              AssemblySystem
        @test AssemblySystem(rules2d, EntropyModel((p1, p2) -> 0.0, COMLaplace())) isa AssemblySystem
    end

    chain = BindingRules([1 1 2 3; 2 1 3 3; 3 1 4 3], UnitSquare)   # a tree: n particles, n-1 bonds
    ring = BindingRules([1 1 2 3; 2 2 3 4; 3 3 4 1; 4 4 1 2], UnitSquare)
    rsp = RigidSpringPotential(0.6; k=120.0)
    ω = bondvolume(rsp)

    # The tree approximation is exactly one bond volume per bond of a tree.
    asys = AssemblySystem(chain, EntropyModel(rsp, TreeLike()); maxsize=4, verbose=false)
    for (s, Ω) in zip(structures(asys), partitionfunctions(asys))
        @test Ω ≈ ω^(nparticles(s) - 1) * 2π / symmetrynumber(s)
        @test length(collect(bondcolors(s))) == nbonds(s)
        @test all(c -> c[1] in 1:ncolors(chain) && c[2] in 1:ncolors(chain), bondcolors(s))
    end
    # and `omega` is superseded once there is a potential to take it from
    @test partitionfunctions(AssemblySystem(chain, EntropyModel(rsp, TreeLike(7.0));
                                            maxsize=4, verbose=false)) ≈ partitionfunctions(asys)
    # a ring still gets charged n-1 bonds, not its n bonds
    asysring = AssemblySystem(ring, EntropyModel(rsp, TreeLike()); maxsize=4, verbose=false)
    Ωlaplace = partitionfunctions(AssemblySystem(ring, EntropyModel(rsp, TetheredLaplace()); maxsize=4, verbose=false))
    @test all(Ω ≈ ω^(nparticles(s) - 1) * 2π / symmetrynumber(s)
              for (s, Ω) in zip(structures(asysring), partitionfunctions(asysring)))

    # The stiff limit of Φ_d is the same Gaussian the Laplace solvers do, so on a tree the two must
    # agree identically rather than just asymptotically.
    for k in (10.0, 1000.0)
        p = RigidSpringPotential(0.7; k)
        tree = partitionfunctions(AssemblySystem(chain, EntropyModel(p, TreeLike());
                                                 maxsize=4, verbose=false))
        laplace = partitionfunctions(AssemblySystem(chain, EntropyModel(p, TetheredLaplace());
                                                    maxsize=4, verbose=false))
        @test tree ≈ laplace rtol = 1e-10
    end
    # …but a closed ring is charged one bond too few, so the tree approximation overshoots there
    let closed = findall(s -> nbonds(s) > nparticles(s) - 1, structures(asysring))
        @test !isempty(closed)
        @test all(partitionfunctions(asysring)[closed] .> Ωlaplace[closed])
        @test partitionfunctions(asysring)[setdiff(1:end, closed)] ≈ Ωlaplace[setdiff(1:end, closed)]
    end

    # `bondvolume` integrates over every relative pose, so a bond that rests elsewhere is fine by
    # it; the Laplace solvers expand about the bonded pose and have to refuse.
    let offmin = RigidSpringPotential(fill([0.0 0.3 -0.2; -0.4 0.25 0.15], ncolors(chain)); k=50.0)
        @test isfinite(bondvolume(offmin))
        @test partitionfunctions(AssemblySystem(chain, EntropyModel(offmin, TreeLike());
                                                maxsize=3, verbose=false)) |> x -> all(isfinite, x)
        for solver in (COMLaplace(), TetheredLaplace())
            @test_throws ArgumentError partitionfunctions(AssemblySystem(chain,
                                                                         EntropyModel(offmin, solver);
                                                                         maxsize=3, verbose=false))
        end
    end

    # Mean field charges each particle for the spread of all of its own patches, so the middle of a
    # chain of squares picks up the half unit between its two sites on top of its clouds.
    let σ = 0.6, k = 120.0, dtot = 3
        mf = AssemblySystem(chain, EntropyModel(RigidSpringPotential(σ; k), MeanField());
                            maxsize=3, verbose=false)
        ω(kk, λ) = sqrt(2π / kk)^dtot / sqrt(λ[1] + λ[2])
        ωend, ωmid = ω(k, [0.0, σ^2 / 12]), ω(2k, [0.25, σ^2 / 12])   # the square's inradius is 1/2
        for (s, Ω) in zip(structures(mf), partitionfunctions(mf))
            n = nparticles(s)
            ωs = n == 1 ? 1.0 : n == 2 ? ωend^2 : ωend^2 * ωmid
            @test Ω ≈ 2π / symmetrynumber(s) * ωs^((n - 1) / n)
        end
    end

    # It is a variational bound, so it can only ever undercount; a dimer it gets exactly right,
    # since tethering one of the two leaves a single free particle and the ansatz costs nothing.
    for rules in (chain, ring)
        mf = partitionfunctions(AssemblySystem(rules, EntropyModel(rsp, MeanField());
                                               maxsize=4, verbose=false))
        exact = AssemblySystem(rules, EntropyModel(rsp, TetheredLaplace()); maxsize=4, verbose=false)
        @test all(mf .<= partitionfunctions(exact) .+ 1e-12)
        for (s, Ωm, Ωl) in zip(structures(exact), mf, partitionfunctions(exact))
            nparticles(s) == 1 && @test Ωm ≈ Ωl                      # nothing to average over
            nparticles(s) == 2 && @test Ωm ≈ Ωl
            nparticles(s) > 2 && @test Ωm < Ωl                       # the ansatz starts to cost
        end
    end

    # ω_i ∝ k_i^(-dtot/2) at every particle, so the structure still scales as k^(-dtot(n-1)/2)
    let a(k) = AssemblySystem(chain, EntropyModel(RigidSpringPotential(0.6; k), MeanField());
                              maxsize=4, verbose=false)
        ns = nparticles.(structures(a(1.0)))
        @test partitionfunctions(a(1.0)) ./ partitionfunctions(a(64.0)) ≈ 64.0 .^ (3 * (ns .- 1) / 2)
    end

    # Mean field is only worked out for a rigid spring potential, on an unstressed structure.
    @test_throws ArgumentError partitionfunctions(AssemblySystem(chain,
                                                                 EntropyModel((p1, p2) -> 0.0,
                                                                              MeanField())))
    let offmin = RigidSpringPotential(fill([0.0 0.3 -0.2; -0.4 0.25 0.15], ncolors(chain)); k=50.0)
        @test_throws ArgumentError partitionfunctions(AssemblySystem(chain,
                                                                     EntropyModel(offmin,
                                                                                  MeanField());
                                                                     maxsize=3, verbose=false))
    end

    # COMLaplace normalises the centre of mass differently, but must still order structures the same
    for solver in (COMLaplace(), TetheredLaplace())
        Ωs = partitionfunctions(AssemblySystem(chain, EntropyModel(rsp, solver); maxsize=4, verbose=false))
        @test all(isfinite, Ωs) && all(>(0), Ωs)
        stiffer = partitionfunctions(AssemblySystem(chain,
                                                    EntropyModel(RigidSpringPotential(0.6; k=480.0),
                                                                 solver); maxsize=4, verbose=false))
        @test all(stiffer .<= Ωs .+ 1e-12)
    end
end;
