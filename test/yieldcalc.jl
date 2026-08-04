begin
    M = [1 0 0;
         0 1 0;
         1 1 1]
    Ωs = ones(3)

    infval = 999

    @test Crafts.nspecies(M) == 2

    @test densities([0, 0, 0], M, Ωs) ≈ Ωs
    @test densities([0, -infval, -infval], M, Ωs) ≈ [1, 0, 0] .* Ωs
    @test densities([-infval, 0, -infval], M, Ωs) ≈ [0, 1, 0] .* Ωs
    @test densities([-infval, -infval, 0], M, Ωs) ≈ [0, 0, 0] .* Ωs
    @test densities([0, 0, -infval], M, Ωs) ≈ [1, 1, 0] .* Ωs

    @test yields([0, 0, 0], M, Ωs) ≈ 1/3 * [1, 1, 1]
    @test yields([0, -infval, -infval], M, Ωs) ≈ [1, 0, 0]
    @test yields([-infval, 0, -infval], M, Ωs) ≈ [0, 1, 0]
    @test yields([0, 0, -infval], M, Ωs) ≈ 1/2 * [1, 1, 0]

    @test particledensities([0, 0, 0], M, Ωs) ≈ [2, 2]
    @test particledensities([0, -infval, -infval], M, Ωs) ≈ [1, 0]
    @test particledensities([-infval, 0, -infval], M, Ωs) ≈ [0, 1]
    @test particledensities([0, 0, -infval], M, Ωs) ≈ [1, 1]

    ϕs = ([1, 0], [0, 1], [1, 1], [0.325, 0.123])
    εs = ([0.0], [1.32])
    for ϕ in ϕs, ε in εs
        μs = chemicalpotentials(ϕ, ε, M, Ωs; atol=1e-8)
        ξ = [μs; only(ε)]
        @test particledensities(ξ, M, Ωs) ≈ ϕ atol=1e-8
    end

    # @test_throws ArgumentError chemical_potentials([0, 0], [1], M, Ωs)
end