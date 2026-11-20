@testset "harmonicpotential" begin
    rng = MersenneTwister(20261119)
    randpose(::Val{2}) = Pose(SVector{2}(randn(rng, 2)), Angle2d(2π * rand(rng)))
    randpose(::Val{3}) = Pose(SVector{3}(randn(rng, 3)), RotXYZ(2π .* rand(rng, 3)...))
    dtot(D) = D * (D + 1) ÷ 2
    randspd(n) = let A = randn(rng, n, n)
        A * A' + n * I
    end
    exchangesymmetric(n) = HarmonicPotential(randspd(n); symmetrize=true).H

    # The coordinates are the ones the Laplace solvers displace a site by, to first order, so the
    # Hessian handed in is the Hessian they see.
    for D in (2, 3)
        J = ForwardDiff.jacobian(ξ -> Crafts._relativecoordinates(Crafts._bondposes(Val(D), ξ)...),
                                 zeros(dtot(D)))
        @test J ≈ I atol = 1e-12
        H = exchangesymmetric(dtot(D))
        @test HarmonicPotential(HarmonicPotential(H); d=D).H ≈ H rtol = 1e-10
    end

    # The energy depends only on the relative pose and not on the order of the two sites.
    for D in (2, 3)
        hp = HarmonicPotential(exchangesymmetric(dtot(D)))
        p1, p2, g = randpose(Val(D)), randpose(Val(D)), randpose(Val(D))
        @test hp(p1, p2) ≈ hp(p2, p1)
        @test hp(g * p1, g * p2) ≈ hp(p1, p2)
        @test hp(p1, Pose(p1.x, p1.psi * Roly.standard_rotation(Float64, Val(D)))) ≈ 0 atol = 1e-12
    end

    # ω is the Gaussian integral of the quadratic form, which the generic expansion reproduces.
    for D in (2, 3)
        H = exchangesymmetric(dtot(D))
        hp = HarmonicPotential(H)
        @test logbondvolume(hp) ≈ dtot(D) / 2 * log(2π) - logdet(H) / 2
        @test logbondvolume((p1, p2, c1, c2) -> hp(p1, p2, c1, c2); d=D) ≈ logbondvolume(hp)
        @test_throws ArgumentError logbondvolume(hp; d=5 - D)
    end

    # The Hessian taken from a function is its second derivative in the relative pose of the sites,
    # at the pose Roly bonds them in.
    let rsp = RigidSpringPotential(0.5, 0.1; k=400.0)
        V(ξ) = rsp(Crafts._bondposes(Val(3), ξ)...)
        @test HarmonicPotential(rsp; d=3).H ≈ ForwardDiff.hessian(V, zeros(6)) rtol = 1e-12
        @test HarmonicPotential((p1, p2, c1, c2) -> rsp(p1, p2); d=3).H ≈ HarmonicPotential(rsp; d=3).H
    end

    # Laplace on a potential and on the harmonic bond built from its own contact Hessian agree for
    # every solver, on trees and on rings alike. The rigid spring has closed forms for `TreeLike` and
    # `MeanField`, so there the closed form is checked against the generic expansion as well.
    chain = BindingRules([1 1 2 3; 2 1 3 3; 3 1 4 3], UnitSquare)
    ring = BindingRules([1 1 2 3; 2 2 3 4; 3 3 4 1; 4 4 1 2], UnitSquare)
    solvers = (COMLaplace(), TetheredLaplace(), MeanField(), TreeLike())
    let rsp = RigidSpringPotential(0.6; k=120.0), hp = HarmonicPotential(rsp; d=2)
        @test logbondvolume(hp) ≈ logbondvolume(rsp)
        for rules in (chain, ring), solver in solvers
            model(p) = AssemblySystem(rules, EntropyModel(p, solver); maxsize=4, verbose=false)
            @test partitionfunctions(model(hp)) ≈ partitionfunctions(model(rsp)) rtol = 1e-8
        end
    end

    # The same holds in 3d, for triangles embedded in space that close into a hexagon.
    hexagon = BindingRules([1 2 2 3; 1 3 2 2], UnitTriangle)
    let rsp = RigidSpringPotential(0.5, 0.1; k=400.0), hp = HarmonicPotential(rsp; d=3)
        model(p, solver) = AssemblySystem(hexagon, EntropyModel(p, solver; embed3d=true); verbose=false)
        @test any(s -> nbonds(s) == nparticles(s), structures(model(hp, COMLaplace())))
        for solver in solvers
            @test partitionfunctions(model(hp, solver)) ≈ partitionfunctions(model(rsp, solver)) rtol = 1e-8
        end
    end

    # An anonymous function outside the rigid spring class, here one that is stiffer normal to the
    # plane and carries anharmonic terms, likewise agrees with the harmonic bond built from it.
    normal(p, q) = (p.psi' * (q.x - p.x))[end]
    for (D, rules, rsp) in ((2, ring, RigidSpringPotential(0.6; k=120.0)),
                            (3, hexagon, RigidSpringPotential(0.5, 0.1; k=400.0)))
        f = (p1, p2, c1, c2) -> rsp(p1, p2) + 300.0 * (normal(p1, p2)^2 + normal(p2, p1)^2) +
                                50.0 * sum(abs2, p2.x - p1.x)^2
        hp = HarmonicPotential(f; d=D)
        @test !(hp.H[1:D, 1:D] ≈ hp.H[1, 1] * I)
        @test logbondvolume(hp) ≈ logbondvolume(f; d=D)
        for solver in solvers
            model(p) = AssemblySystem(rules, EntropyModel(p, solver; embed3d=D == 3); maxsize=6, verbose=false)
            @test partitionfunctions(model(hp)) ≈ partitionfunctions(model(f)) rtol = 1e-8
        end
    end

    # A function that is not symmetric under exchanging the sites gives an asymmetric Hessian, which
    # has to be symmetrized explicitly. In 2d exchanging the sites keeps Δx and flips φ, so a coupling
    # between the two is odd.
    let rsp = RigidSpringPotential(0.6; k=120.0),
        skewed = (p1, p2, c1, c2) -> rsp(p1, p2) + 30.0 * prod(Crafts._relativecoordinates(p1, p2)[2:3])
        @test_throws ArgumentError HarmonicPotential(skewed; d=2)
        @test HarmonicPotential(skewed; d=2, symmetrize=true) isa HarmonicPotential
    end

    # A Hessian outside the rigid spring class still gives a finite, well ordered answer, and a
    # stiffer bond can only reduce every partition function.
    let H = exchangesymmetric(6),
        model(H) = AssemblySystem(hexagon, EntropyModel(HarmonicPotential(H), COMLaplace(); embed3d=true);
                                  verbose=false)
        Ω = partitionfunctions(model(H))
        @test all(isfinite, Ω) && all(>(0), Ω)
        @test all(partitionfunctions(model(4H)) .<= Ω)
    end

    # Malformed Hessians are refused.
    @test_throws ArgumentError HarmonicPotential(randspd(4))
    @test_throws ArgumentError HarmonicPotential(randspd(6) + triu(ones(6, 6), 1))
    @test_throws ArgumentError HarmonicPotential(-exchangesymmetric(6))

    # A Hessian that changes when the sites are exchanged is refused unless symmetrized, and
    # symmetrizing is a projection.
    for D in (2, 3)
        H = randspd(dtot(D))
        J = Crafts._exchangejacobian(Val(D), Float64)
        @test J' * J ≈ I
        @test !(H ≈ J' * H * J)
        @test_throws ArgumentError HarmonicPotential(H)
        Hs = HarmonicPotential(H; symmetrize=true).H
        @test Hs ≈ (H + J' * H * J) / 2
        @test HarmonicPotential(Hs).H ≈ Hs
    end

    # A potential of the wrong dimension is caught when the system is built.
    @test_throws ArgumentError AssemblySystem(chain, EntropyModel(HarmonicPotential(exchangesymmetric(6)),
                                                                  TreeLike()))
    @test_throws ArgumentError AssemblySystem(hexagon,
                                              EntropyModel(HarmonicPotential(exchangesymmetric(3)), TreeLike();
                                                           embed3d=true))
end;
