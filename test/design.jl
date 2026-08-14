@testset "design" begin
    rules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false)
    M, Ωs = compositionmatrix(asys), partitionfunctions(asys)
    trimer = 6
    opt = Clarabel.Optimizer

    @testset "convexdesign" begin
        ξ, residual = convexdesign(asys, trimer; maxε=8.0, optimizer=opt)
        ys = yields(ξ, M, Ωs)
        # designing at the same bond energy must beat handing out that energy uniformly
        @test ys[trimer] > yields(fill(0.1, 3), fill(8.0, 2), M, Ωs)[trimer]
        @test residual ≈ log(sum(ys[1:(end - 1)]) / ys[trimer]) rtol = 1e-4

        # a scalar target and a one-element vector are the same problem
        ξv, _ = convexdesign(asys, [trimer]; maxε=8.0, optimizer=opt)
        @test ξ ≈ ξv atol = 1e-5
    end

    @testset "relative yields" begin
        ξ, _ = convexdesign(asys, [4, 5]; maxε=8.0, optimizer=opt)
        ρ = densities(ξ, M, Ωs)
        @test ρ[4] ≈ ρ[5] rtol = 1e-4

        ξ2, _ = convexdesign(asys, [4, 5]; relative_yields=[2.0, 1.0], maxε=8.0, optimizer=opt)
        ρ2 = densities(ξ2, M, Ωs)
        @test ρ2[4] / ρ2[5] ≈ 2 rtol = 1e-4
    end

    @testset "minenergydesign inverts convexdesign" begin
        for y in (0.5, 0.7, 0.9)
            ξ, ε̄ = minenergydesign(asys, trimer; minyield=y, optimizer=opt)
            @test yields(ξ, M, Ωs)[trimer] ≈ y rtol = 1e-3

            # feeding the achieved energy back must recover the same yield
            ξback, _ = convexdesign(asys, trimer; maxε=ε̄, optimizer=opt)
            @test yields(ξback, M, Ωs)[trimer] ≈ y rtol = 1e-3
        end
    end

    @testset "εbound" begin
        ξmean, _ = minenergydesign(asys, trimer; minyield=0.8, εbound=:mean, optimizer=opt)
        ξmax, _ = minenergydesign(asys, trimer; minyield=0.8, εbound=:max, optimizer=opt)
        # minimizing the largest bond energy cannot leave a larger one than minimizing the mean
        @test maximum(ξmax[4:5]) <= maximum(ξmean[4:5]) + 1e-4
        @test yields(ξmax, M, Ωs)[trimer] ≈ 0.8 rtol = 1e-3
    end

    @testset "minenergydesign with several targets" begin
        # the yield floor applies to the target set as a whole, so the bound on the competition
        # scales with how much of the total the targets already claim
        for (idxs, ry, y) in (([4, 5], [1.0, 1.0], 0.4), ([4, 5], [3.0, 1.0], 0.4),
                              ([4, 5, 6], [1.0, 1.0, 2.0], 0.75))
            ξ, _ = minenergydesign(asys, idxs; minyield=y, relative_yields=ry, optimizer=opt)
            ys = yields(ξ, M, Ωs)
            @test sum(ys[idxs]) ≈ y rtol = 1e-3
            @test ys[idxs] ./ ys[idxs[1]] ≈ ry ./ ry[1] rtol = 1e-3
        end
    end

    @testset "infeasible yield" begin
        @test_throws ArgumentError minenergydesign(asys, trimer; minyield=1.0, optimizer=opt)
        @test_throws ArgumentError minenergydesign(asys, trimer; minyield=0.0, optimizer=opt)

        # the two dimers compete with each other and with the trimer, so they cannot hold most of
        # the population between them; that must be reported rather than silently returning a point
        # that violates the constraints
        @test_throws ArgumentError minenergydesign(asys, [4, 5]; minyield=0.9, optimizer=opt)
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
