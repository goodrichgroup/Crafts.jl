@testset "design" begin
    rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false)
    M, Ωs = compositionmatrix(asys), partitionfunctions(asys)
    trimer = 6
    opt = Clarabel.Optimizer

    @testset "maxyielddesign" begin
        ξ, residual = maxyielddesign(asys, trimer; maxdensity=1, energy_budget=8.0, optimizer=opt)
        ys = yields(ξ, M, Ωs)
        # designing at the same bond energy must beat handing out that energy uniformly
        @test ys[trimer] > yields(fill(0.1, 3), fill(8.0, 2), M, Ωs)[trimer]
        @test residual ≈ log(sum(ys[1:(end - 1)]) / ys[trimer]) rtol = 1e-4
        # the budget is a cap on the mean, and here it pays to spend all of it
        @test mean(ξ[4:5]) ≈ 8.0 rtol = 1e-4

        # a scalar target and a one-element vector are the same problem
        ξv, _ = maxyielddesign(asys, [trimer]; maxdensity=1, energy_budget=8.0, optimizer=opt)
        @test ξ ≈ ξv atol = 1e-5

        # the density scale sets what the yields mean, so it is stated rather than defaulted
        @test_throws UndefKeywordError maxyielddesign(asys, trimer; energy_budget=8.0, optimizer=opt)
        @test_throws UndefKeywordError minenergydesign(asys, trimer; minyield=0.5, optimizer=opt)
    end

    @testset "relative yields" begin
        ξ, _ = maxyielddesign(asys, [4, 5]; maxdensity=1, energy_budget=8.0, optimizer=opt)
        ρ = densities(ξ, M, Ωs)
        @test ρ[4] ≈ ρ[5] rtol = 1e-4

        ξ2, _ = maxyielddesign(asys, [4, 5]; maxdensity=1, relative_yields=[2.0, 1.0], energy_budget=8.0,
                               optimizer=opt)
        ρ2 = densities(ξ2, M, Ωs)
        @test ρ2[4] / ρ2[5] ≈ 2 rtol = 1e-4
    end

    @testset "minenergydesign inverts maxyielddesign" begin
        for y in (0.5, 0.7, 0.9)
            ξ, ε̄ = minenergydesign(asys, trimer; maxdensity=1, minyield=y, optimizer=opt)
            @test yields(ξ, M, Ωs)[trimer] ≈ y rtol = 1e-3

            # feeding the achieved energy back must recover the same yield
            ξback, _ = maxyielddesign(asys, trimer; maxdensity=1, energy_budget=ε̄, optimizer=opt)
            @test yields(ξback, M, Ωs)[trimer] ≈ y rtol = 1e-3
        end
    end

    @testset "energy_measure" begin
        ξmean, _ = minenergydesign(asys, trimer; maxdensity=1, minyield=0.8, energy_measure=mean, optimizer=opt)
        ξmax, r = minenergydesign(asys, trimer; maxdensity=1, minyield=0.8, energy_measure=maximum, optimizer=opt)
        # minimizing the largest bond energy cannot leave a larger one than minimizing the mean
        @test maximum(ξmax[4:5]) <= maximum(ξmean[4:5]) + 1e-4
        @test yields(ξmax, M, Ωs)[trimer] ≈ 0.8 rtol = 1e-3
        # the residual is the measure that was minimized, whichever one it is
        @test r ≈ maximum(ξmax[4:5]) rtol = 1e-4

        # the budget bounds whatever the measure reports, so capping the largest bond energy is the
        # tighter request: it cannot spend more in the mean than a mean budget would allow
        ξb, _ = maxyielddesign(asys, trimer; maxdensity=1, energy_budget=8.0, energy_measure=maximum,
                               optimizer=opt)
        @test maximum(ξb[4:5]) <= 8.0 + 1e-4
    end

    @testset "uniform_energy" begin
        ξ, _ = maxyielddesign(asys, trimer; maxdensity=1, energy_budget=8.0, uniform_energy=true, optimizer=opt)
        @test ξ[4] ≈ ξ[5] rtol = 1e-5
        # one energy for every bond means mean and maximum coincide, and both sit at the budget
        @test ξ[4] ≈ 8.0 rtol = 1e-4

        ξu, εu = minenergydesign(asys, trimer; maxdensity=1, minyield=0.8, uniform_energy=true, optimizer=opt)
        ξf, εf = minenergydesign(asys, trimer; maxdensity=1, minyield=0.8, optimizer=opt)
        @test ξu[4] ≈ ξu[5] rtol = 1e-5
        @test yields(ξu, M, Ωs)[trimer] ≈ 0.8 rtol = 1e-3
        # tying the bonds together can only cost energy, never save it
        @test εu >= εf - 1e-4
    end

    @testset "uniform_energy over three bond types" begin
        # two bond types tie with a single equality; three make sure the whole chain is tied
        rules3 = BindingRules([1 1 2 3; 2 1 3 3; 3 1 4 3], UnitSquare)
        asys3 = AssemblySystem(rules3, EntropyModel(TreeLike()); verbose=false)
        M3, Ω3 = compositionmatrix(asys3), partitionfunctions(asys3)
        chain = nstructures(asys3)                          # the full four-particle chain
        bonds = 5:7

        ξ, _ = maxyielddesign(asys3, chain; maxdensity=1, energy_budget=6.0, uniform_energy=true, optimizer=opt)
        @test all(≈(ξ[first(bonds)]; rtol=1e-5), ξ[bonds])
        @test ξ[first(bonds)] ≈ 6.0 rtol = 1e-4

        # left free, the middle bond takes more than the outer two, so the tie is not vacuous
        ξf, _ = maxyielddesign(asys3, chain; maxdensity=1, energy_budget=6.0, optimizer=opt)
        @test !isapprox(ξf[6], ξf[5]; rtol=1e-3)
        # a mean budget bounds only the average, so an unequal solution must exceed it somewhere
        @test mean(ξf[bonds]) ≈ 6.0 rtol = 1e-4
        @test maximum(ξf[bonds]) > 6.0

        ξu, _ = minenergydesign(asys3, chain; maxdensity=1, minyield=0.5, uniform_energy=true, optimizer=opt)
        @test all(≈(ξu[first(bonds)]; rtol=1e-5), ξu[bonds])
        @test yields(ξu, M3, Ω3)[chain] ≈ 0.5 rtol = 1e-3
    end

    @testset "uniform_energy with one bond type" begin
        # nothing to tie, so the constraint has to be skipped rather than built over an empty range
        rules1 = BindingRules([1 1 2 1], UnitSquare)
        asys1 = AssemblySystem(rules1, EntropyModel(TreeLike()); verbose=false)
        dimer = nstructures(asys1)

        ξ, residual = maxyielddesign(asys1, dimer; maxdensity=1, energy_budget=5.0, uniform_energy=true,
                                     optimizer=opt)
        ξf, residualf = maxyielddesign(asys1, dimer; maxdensity=1, energy_budget=5.0, optimizer=opt)
        @test ξ ≈ ξf atol = 1e-6
        @test residual ≈ residualf rtol = 1e-6
    end

    @testset "minenergydesign with several targets" begin
        # the yield floor applies to the target set as a whole, so the bound on the competition
        # scales with how much of the total the targets already claim
        for (idxs, ry, y) in (([4, 5], [1.0, 1.0], 0.4), ([4, 5], [3.0, 1.0], 0.4),
                              ([4, 5, 6], [1.0, 1.0, 2.0], 0.75))
            ξ, _ = minenergydesign(asys, idxs; maxdensity=1, minyield=y, relative_yields=ry, optimizer=opt)
            ys = yields(ξ, M, Ωs)
            @test sum(ys[idxs]) ≈ y rtol = 1e-3
            @test ys[idxs] ./ ys[idxs[1]] ≈ ry ./ ry[1] rtol = 1e-3
        end
    end

    @testset "infeasible yield" begin
        @test_throws ArgumentError minenergydesign(asys, trimer; maxdensity=1, minyield=1.0, optimizer=opt)
        @test_throws ArgumentError minenergydesign(asys, trimer; maxdensity=1, minyield=0.0, optimizer=opt)

        # the two dimers compete with each other and with the trimer, so they cannot hold most of
        # the population between them; that must be reported rather than silently returning a point
        # that violates the constraints
        @test_throws ArgumentError minenergydesign(asys, [4, 5]; maxdensity=1, minyield=0.9, optimizer=opt)
    end

    @testset "lineardesign" begin
        ξ, _ = lineardesign(asys, trimer; optimizer=opt)
        @test length(ξ) == size(M, 2)
        # the target must sit at zero excess free energy, and nothing else may join it
        @test isapprox(dot(M[trimer, :], ξ), 0; atol=1e-5)
        @test all(<(0), M[setdiff(1:end, trimer), :] * ξ .- 1e-5)
    end

    @testset "preprocessing" begin
        mask, elements, newidxs = Crafts._preprocessdesign(M, [4])
        # the dimer of species 1 and 2 uses neither species 3 nor the second bond
        @test !elements[3] && !elements[5]
        @test mask == [true, true, false, true, false, false]
        @test newidxs == [3]
    end
end
