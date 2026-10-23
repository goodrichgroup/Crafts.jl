@testset "bindingruledesign" begin
    ### A 2x2 block of four distinct squares, sites 1,2,3,4 running around each square:
    ### particle 1 top-left, 2 top-right, 3 bottom-right, 4 bottom-left.
    sq = PolygonParticleSpecies(4, 1.0)
    rules = BindingRules([1 2 2 4; 2 3 3 1; 3 4 4 2; 4 1 1 3], sq)

    strs = polygen(rules; maxsize=8, maxpolys=500)
    target = only(filter(s -> nparticles(s) == 4, strs))
    @test nbonds(target) == 4

    canonical(r) = first(Crafts._canonicalforms(r))

    @testset "fully_addressable_design" begin
        fad = fully_addressable_design(target, sq)
        @test nspecies(fad) == 4
        @test nbonds(fad) == 4

        # The block was built from four mutually distinct species, so it is already fully
        # addressable and the seed must recover the rules it came from.
        @test canonical(fad) == canonical(rules)
        @test length(polygen(fad; maxsize=8, maxpolys=500)) == length(strs)

        @test_throws ArgumentError fully_addressable_design(target, fill(sq, 3))
    end

    @testset "canonical forms" begin
        fad = fully_addressable_design(target, sq)
        form, mirror = Crafts._canonicalforms(fad)
        # Both handednesses are stored so that one lookup recognizes either.
        @test form in Set(Crafts._canonicalforms(fad))
        @test mirror in Set(Crafts._canonicalforms(fad))

        # Merging changes the rules, so the canonical form has to change with it.
        @test canonical(fad) != canonical(Crafts.merge_species(fad, 1, 2, 0))
    end

    @testset "find_shape" begin
        @test find_shape(target, strs) == [length(strs)]
        # Shape alone, so the monomers all match each other regardless of species.
        monomer = first(strs)
        @test length(find_shape(monomer, strs)) == count(s -> nparticles(s) == 1, strs)
    end

    @testset "merge_species and refinements" begin
        fad = fully_addressable_design(target, sq)
        merged = Crafts.merge_species(fad, 1, 2, 0)
        @test nspecies(merged) == nspecies(fad) - 1

        # Merging only ever adds bonds, so no structure is lost.
        @test !isempty(find_shape(target, polygen(merged; maxsize=8, maxpolys=500)))

        # Six unordered pairs of the four species, four relative rotations each.
        @test length(Crafts.generate_refinements(fad)) == 24
        @test length(Crafts.generate_refinements(fad; enable_rotations=false)) == 6
    end

    @testset "evaluate_system" begin
        fad = fully_addressable_design(target, sq)
        valid, score = Crafts.evaluate_system(fad; maxsize=8, maxpolys=500, target)
        @test valid
        @test score == -length(strs)

        # Too tight a bound makes the rule set look unbounded, so it is rejected.
        @test first(Crafts.evaluate_system(fad; maxsize=3, maxpolys=500, target)) == false
        @test first(Crafts.evaluate_system(fad; maxsize=8, maxpolys=5, target)) == false

        valid, yield = Crafts.evaluate_system(fad, target, -12.0, 1e-3; maxsize=8, maxpolys=500)
        @test valid
        @test 0 <= yield <= 1
    end

    @testset "create_designs" begin
        systems, network = create_designs(target, sq; maxsize=8, maxpolys=500)

        @test length(systems) == Graphs.nv(network)
        # The search starts at the seed and keeps it as the first entry.
        @test canonical(first(systems)) == canonical(fully_addressable_design(target, sq))

        # Every design must still assemble the target, and none may repeat.
        for sys in systems
            @test !isempty(find_shape(target, polygen(sys; maxsize=8, maxpolys=500)))
        end
        forms = canonical.(systems)
        @test length(Set(forms)) == length(forms)

        # Merging strictly reduces the species count, so the search reaches a single species.
        @test minimum(nspecies, systems) == 1

        yieldscored, _ = create_designs(target, sq; maxsize=8, maxpolys=500, ε=-12.0, ϕ=1e-3)
        @test !isempty(yieldscored)
        for sys in yieldscored
            @test !isempty(find_shape(target, polygen(sys; maxsize=8, maxpolys=500)))
        end
    end

    @testset "create_multifarious_designs" begin
        fad = fully_addressable_design(target, sq)
        systems, _ = create_multifarious_designs(fad, [target]; maxsize=8, maxpolys=500)
        @test !isempty(systems)
        for sys in systems
            @test !isempty(find_shape(target, polygen(sys; maxsize=8, maxpolys=500)))
        end
    end
end
