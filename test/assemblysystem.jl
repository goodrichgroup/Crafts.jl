@testset "assemblysystem" begin
    rules = BindingRules(
        [1 3 2 3;
         2 1 4 1;
         2 2 3 2;
         3 1 4 1],
        UnitTriangle)   # 16 structures, largest has 5 particles
    model = EntropyModel(TreeLike())

    asys = AssemblySystem(rules, model)
    @test nspecies(asys) == 4
    @test nbonds(asys) == 4
    @test dimension(asys) == 2
    @test asys.maxsize == Inf
    @test ismissing(iscomplete(asys))

    # counting must not retain the structures themselves
    @test nstructures(asys) == 16
    @test asys.M !== nothing
    @test asys.structures === nothing
    @test size(compositionmatrix(asys)) == (16, nspecies(asys) + nbonds(asys))
    @test iscomplete(asys)

    strs = structures(asys)
    @test length(strs) == 16
    @test issorted(nparticles.(strs))
    # `M` must be built from this exact ordering, or it would misalign with `Ωs`
    @test compositionmatrix(asys) == reduce(vcat, transpose.(composition.(strs)))

    Ωs = partitionfunctions(asys)
    @test length(Ωs) == 16
    @test Ωs ≈ [2π / symmetrynumber(p) for p in strs]        # omega = 1

    half = AssemblySystem(rules, EntropyModel(TreeLike(0.5)))
    @test partitionfunctions(half) ≈ [0.5^(nparticles(p) - 1) * 2π / symmetrynumber(p) for p in structures(half)]

    # solvers that need a potential cannot be built without one
    @test_throws ArgumentError EntropyModel(COMLaplace())
    @test_throws ArgumentError EntropyModel(TetheredLaplace())

    # check warning for incomplete enumerations
    let cut = AssemblySystem(rules, model; maxsize=3)
        @test (@test_logs (:warn,) match_mode = :any structures(cut)) === structures(cut)
        @test iscomplete(cut) == false
        @test_logs structures(cut)               # cached, so silent
        @test_logs compositionmatrix(cut)
        @test_logs partitionfunctions(cut)
    end
    let whole = AssemblySystem(rules, model)
        @test_logs structures(whole)
        @test iscomplete(whole)
    end

    # a narrower one-off must not disturb the cache
    @test nstructures(asys; maxsize=3) == 12
    @test size(asys.M, 1) == 16
    @test length(structures(asys; maxsize=3)) == 12
    @test length(asys.structures) == 16
    # asking for the system's own maxsize explicitly still hits the cache
    @test nstructures(asys; maxsize=Inf) == 16

    capped = AssemblySystem(rules, model; maxsize=3, verbose=false)
    @test capped.maxsize == 3
    @test nstructures(capped) == 12
    @test all(≤(3), nparticles.(structures(capped)))
    @test !iscomplete(capped)

    # a cutoff above the largest structure does not actually cut anything
    loose = AssemblySystem(rules, model; maxsize=100, verbose=false)
    @test nstructures(loose) == 16
    @test iscomplete(loose)

    @test_throws ArgumentError AssemblySystem(rules, model; maxsize=0)
    @test_throws ArgumentError AssemblySystem(rules, model; maxsize=-1)
    @test_throws ArgumentError AssemblySystem(rules, model; maxsize=2.5)

    c = countstructures(asys)
    @test c.exact
    @test c.n == 16
    @test c.largest_size == 5
    @test countstructures(asys; maxsize=3).n == 12

    @test occursin("n=4, k=4", sprint(show, asys))
    @test occursin("maxsize=3", sprint(show, capped))
end;
