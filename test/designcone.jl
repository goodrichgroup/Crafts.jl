@testset "designcone" begin
    # The minimal example of Hübl et al., Nat. Phys. (2026), Fig. 3a: a threefold-symmetric central
    # triangle capped by up to three copies of a second species. Its three sites share a colour,
    # which is what makes them equivalent, and so they are one bond type carrying one energy.
    central = PolygonParticleSpecies(3; colors=[1, 1, 1])
    outer = PolygonParticleSpecies(3)
    rules = BindingRules([1 1 2 1], [central, outer])
    asys = AssemblySystem(rules, EntropyModel(TreeLike()); verbose=false)
    CENTRAL, OUTER, DIMER, TRIMER, TETRAMER = 1, 2, 3, 4, 5

    # The same capping geometry with the central sites left distinguishable, so it carries three
    # bond types instead of one. Grouping has nothing to do on a single bond type, so the tests
    # about grouping use this rather than the figure's system.
    asym = AssemblySystem(BindingRules([1 1 2 1; 1 2 2 1; 1 3 2 1],
                                       [PolygonParticleSpecies(3), PolygonParticleSpecies(3)]),
                          EntropyModel(TreeLike()); verbose=false)

    @testset "parametermap" begin
        # the figure's system has one bond type, so its map is already the identity
        @test size(parametermap(asys)) == (3, 3)
        @test parametermap(asys) == I
        @test parametermap(asys; bondgroups=:uniform) == I

        # three bond types, where grouping does something
        @test size(parametermap(asym)) == (5, 5)
        @test parametermap(asym) == I
        @test size(parametermap(asym; bondgroups=:uniform)) == (5, 3)
        @test Crafts._reducedM(asym; bondgroups=:uniform) ==
              [1 0 0; 0 1 0; 1 1 1; 1 1 1; 1 1 1; 1 2 2; 1 2 2; 1 2 2; 1 3 3]
        @test parametermap(asym; bondgroups=[[1, 2], [3]]) ==
              parametermap(asym; bondgroups=[[2, 1], [3]])
        @test_throws ArgumentError parametermap(asym; bondgroups=[[1, 2]])
        @test_throws ArgumentError parametermap(asym; bondgroups=[[1, 1], [2], [3]])
    end

    c = constraintcone(asys; bondgroups=:uniform)

    @testset "Fig. 3a designability" begin
        @test size(c.M, 2) == 3          # two chemical potentials and one binding energy
        @test size(c.rays, 1) == 3

        # both monomers and the tetramer bound the cone; the dimer and trimer only touch it
        @test designablestructures(asys; bondgroups=:uniform) == [CENTRAL, OUTER, TETRAMER]
        @test isdesignable(c, TETRAMER)
        @test !isdesignable(c, DIMER)
        @test !isdesignable(c, TRIMER)

        # the dimer and trimer can only assemble alongside the tetramer and the central monomer
        @test minimaldesignableset(c, DIMER) == [CENTRAL, DIMER, TRIMER, TETRAMER]
        @test minimaldesignableset(c, TRIMER) == [CENTRAL, DIMER, TRIMER, TETRAMER]
        @test necessarychimeras(c, DIMER) == [CENTRAL, TRIMER, TETRAMER]
        @test isempty(necessarychimeras(c, TETRAMER))

        # the three rays of the cone, as in the paper's Hasse diagram
        sets = designablesets(c)
        @test [CENTRAL, OUTER] in sets
        @test [OUTER, TETRAMER] in sets
        @test [CENTRAL, DIMER, TRIMER, TETRAMER] in sets
    end

    @testset "tying bond energies restricts design" begin
        # with every bond independent there is enough freedom to isolate any structure
        @test designablestructures(asym) == 1:9
        @test isdesignable(constraintcone(asym), 3)
        # tie the three energies together and only the monomers and the full tetramer survive
        @test designablestructures(asym; bondgroups=:uniform) == [1, 2, 9]
    end

    @testset "designable set invariants" begin
        for cone in (c, constraintcone(asym))
            sets = designablesets(cone)
            n = size(cone.M, 1)

            @test sort(vcat([s for s in sets if length(s) == 1]...)) == sort(cone.designable)
            @test all(s -> isempty(s) || isdesignable(cone, s), sets)
            @test length(unique(sort.(sets))) == length(sets)

            Random.seed!(1)
            for _ in 1:50
                T = sort(randperm(n)[1:rand(1:min(4, n))])
                S = minimaldesignableset(cone, T)
                @test issubset(T, S)
                @test minimaldesignableset(cone, S) == S   # closure is idempotent
            end

            # designable sets are closed under intersection, so a minimal one always exists
            for a in sets, b in sets
                x = intersect(a, b)
                isempty(x) || @test isdesignable(cone, x)
            end
        end
    end

    @testset "maxsize truncation" begin
        full = designablesets(c)
        for k in 1:4
            @test designablesets(c; maxsize=k) == filter(s -> length(s) <= k, full)
        end
    end

    @testset "reduced design space" begin
        set = [CENTRAL, DIMER, TRIMER, TETRAMER]
        @test codimension(asys, set; bondgroups=:uniform) == 2
        @test codimension(c, set) == 2
        # one of the two directions only rescales the whole set, leaving one to tune yields
        @test relativeyielddofs(asys, set; bondgroups=:uniform) == 1
        @test size(yielddirections(asys, set; bondgroups=:uniform)) == (3, 1)

        # at the origin the relative yields are the bare partition functions, which for three
        # equivalent binding sites are the binomial multiplicities
        y = relativeyields(asys, set, [0.0]; bondgroups=:uniform)
        @test y ≈ [1, 3, 3, 1] ./ 8 rtol = 1e-8
        @test sum(y) ≈ 1

        # sweeping the free direction moves the population from the monomer to the tetramer
        @test relativeyields(asys, set, [-6.0]; bondgroups=:uniform)[1] > 0.99
        @test relativeyields(asys, set, [6.0]; bondgroups=:uniform)[end] > 0.99

        @test_throws DimensionMismatch relativeyields(asys, set, [0.0, 0.0]; bondgroups=:uniform)
    end

    @testset "uniform direction does not change relative yields" begin
        set = [CENTRAL, DIMER, TRIMER, TETRAMER]
        MT = Float64.(Crafts._reducedM(asys; bondgroups=:uniform)[set, :])
        Ω = partitionfunctions(asys)[set]
        u = pinv(MT) * ones(length(set))
        norm(u) > 0 || error("expected a uniform rescaling direction")
        rel(ξ) = (v = Ω .* exp.(MT * ξ); v ./ sum(v))
        @test rel(zeros(3)) ≈ rel(3u) rtol = 1e-10
    end
end
