@testset "yields" begin
    M = [1 0 0;
         0 1 0;
         1 1 1]
    Ωs = ones(3)

    infval = 999

    @test Crafts._nspecies(M) == 2

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
        μsol = chemicalpotentials(ϕ, ε, M, Ωs; atol=1e-8)
        ξsol = [μsol; only(ε)]
        @test particledensities(ξsol, M, Ωs) ≈ ϕ atol=1e-8
    end

    @test_throws ArgumentError chemicalpotentials([0, 0], [1], M, Ωs)

    # two species, one bond type: reproduce the same `M` from above
    rules = BindingRules([1 1 2 1], UnitTriangle)
    asys = AssemblySystem(rules, EntropyModel(TreeApproximation()))

    M = compositionmatrix(asys)
    Ωs = partitionfunctions(asys)
    @test M == [1 0 0; 0 1 0; 1 1 1]

    ξ = [-1.0, -2.0, 3.0]
    @test logdensities(asys, ξ) ≈ logdensities(ξ, M, Ωs)
    @test densities(asys, ξ) ≈ densities(ξ, M, Ωs)
    @test logyields(asys, ξ) ≈ logyields(ξ, M, Ωs)
    @test yields(asys, ξ) ≈ yields(ξ, M, Ωs)
    @test logparticledensities(asys, ξ) ≈ logparticledensities(ξ, M, Ωs)
    @test particledensities(asys, ξ) ≈ particledensities(ξ, M, Ωs)

    @test sum(yields(asys, ξ)) ≈ 1
    @test densities(asys, [0, 0, 0]) ≈ Ωs

    ϕs, εs = [0.3, 0.2], [3.0]
    μs = chemicalpotentials(asys, ϕs, εs)
    @test μs ≈ chemicalpotentials(ϕs, εs, M, Ωs)
    @test particledensities([μs; εs], M, Ωs) ≈ ϕs atol=1e-6
    @test yields(asys, ϕs, εs) ≈ yields(ϕs, εs, M, Ωs)
    @test densities(asys, ϕs, εs) ≈ densities(ϕs, εs, M, Ωs)
    @test logdensities(asys, ϕs, εs) ≈ logdensities(ϕs, εs, M, Ωs)
    @test logyields(asys, ϕs, εs) ≈ logyields(ϕs, εs, M, Ωs)
    @test logparticledensities(asys, ϕs, εs) ≈ logparticledensities(ϕs, εs, M, Ωs)
    @test particledensities(asys, ϕs, εs) ≈ particledensities(ϕs, εs, M, Ωs)

    # every log form is the log of its counterpart
    @test logdensities(asys, ξ) ≈ log.(densities(asys, ξ))
    @test logyields(asys, ξ) ≈ log.(yields(asys, ξ))
    @test logparticledensities(asys, ξ) ≈ log.(particledensities(asys, ξ))
    @test logdensities(asys, ϕs, εs) ≈ log.(densities(asys, ϕs, εs))
    @test logyields(asys, ϕs, εs) ≈ log.(yields(asys, ϕs, εs))
    @test logparticledensities(asys, ϕs, εs) ≈ log.(particledensities(asys, ϕs, εs))

    # check if `solve_kwargs` reach the solver
    @test chemicalpotentials(asys, [0.3, 0.0], εs; infval=42)[2] == -42
    @test chemicalpotentials(asys, [0.3, 0.0], εs)[2] != -42

    # check output eltype
    @test eltype(Ωs) == Float64
    @test eltype(chemicalpotentials(asys, Float32.(ϕs), Float32.(εs))) == Float32
    @test eltype(chemicalpotentials(asys, ϕs, εs)) == Float64
    @test eltype(chemicalpotentials(asys, [1, 0], [3])) == Float64  # integers promote to float

    # check solver algorithm pass through
    @test chemicalpotentials(asys, ϕs, εs; alg=NewtonRaphson()) ≈ chemicalpotentials(asys, ϕs, εs)
    @test chemicalpotentials(asys, ϕs, εs; alg=nothing) ≈ chemicalpotentials(asys, ϕs, εs)

    # other pass through checks
    for f in (logdensities, densities, logyields, yields, logparticledensities, particledensities,
              chemicalpotentials, density_jacobian, yield_jacobian)
        @test f(asys, ϕs, εs; atol=1e-10, rtol=1e-10) ≈ f(ϕs, εs, M, Ωs; atol=1e-10, rtol=1e-10)
        @test_throws MethodError f(asys, ϕs, εs; not_a_solver_option=1)
    end

    @test_throws ArgumentError yields(asys, [0.0, 0.0])
    @test_throws ArgumentError yields(asys, [1.0], εs)
    @test_throws ArgumentError chemicalpotentials(asys, [1.0], εs)
    @test_throws ArgumentError logdensities(asys, [0.0, 0.0])
end;

@testset "jacobians" begin
    rules = BindingRules(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangle)
    asys = AssemblySystem(rules, EntropyModel(TreeApproximation()))
    M = compositionmatrix(asys)
    Ωs = partitionfunctions(asys)
    ns, nb = nspecies(asys), nbonds(asys)

    ξ = [-3.0, -3.5, -4.0, -2.5, 1.0, 1.5, 0.5, 2.0]

    @test density_jacobian(asys, ξ) ≈ ForwardDiff.jacobian(x -> densities(asys, x), ξ)
    @test yield_jacobian(asys, ξ) ≈ ForwardDiff.jacobian(x -> yields(asys, x), ξ)
    @test particledensity_jacobian(asys, ξ) ≈ ForwardDiff.jacobian(x -> particledensities(asys, x), ξ)

    @test density_jacobian(asys, ξ) ≈ density_jacobian(ξ, M, Ωs)
    @test yield_jacobian(asys, ξ) ≈ yield_jacobian(ξ, M, Ωs)
    @test particledensity_jacobian(asys, ξ) ≈ particledensity_jacobian(ξ, M, Ωs; ns)

    # yields are normalised, so their derivatives cancel across structures
    @test all(≈(0; atol=1e-12), sum(yield_jacobian(asys, ξ); dims=1))

    # Check the ϕ-space jacobians against relations they must satisfy
    ϕs, εs = [0.20, 0.15, 0.10, 0.25], [1.0, 1.5, 0.5, 2.0]
    N = view(M, :, 1:ns)
    Jρ = density_jacobian(asys, ϕs, εs)
    @test size(Jρ) == (nstructures(asys), ns + nb)
    @test N' * Jρ[:, 1:ns] ≈ I(ns) atol = 1e-6
    @test N' * Jρ[:, ns+1:end] ≈ zeros(ns, nb) atol = 1e-6
    @test Jρ ≈ density_jacobian(ϕs, εs, M, Ωs)

    JY = yield_jacobian(asys, ϕs, εs)
    @test all(≈(0; atol=1e-10), sum(JY; dims=1))
    @test JY ≈ yield_jacobian(ϕs, εs, M, Ωs)
end;
