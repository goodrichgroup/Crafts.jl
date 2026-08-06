@testset "reactions" begin
    ncuts(p; maxbonds) = length(first(Crafts.generate_cuts(graphrep(p); maxbonds)))
    target(rules, maxsize) = argmax(nbonds,
        structures(AssemblySystem(rules, EntropyModel(TreeApproximation()); maxsize, verbose=false)))

    # Rings of length 4, 6, 6 and 8, built from squares, triangles, and hexagons
    rings = ((BindingRules([1 1 2 3; 2 2 3 4; 3 3 4 1; 4 4 1 2], UnitSquare), 4),
             (BindingRules([1 1 2 1; 2 2 3 2; 3 3 4 3; 4 1 5 1; 5 2 6 2; 6 3 1 3], UnitTriangle), 6),
             (BindingRules([1 1 2 4; 2 2 3 5; 3 3 4 6; 4 4 5 1; 5 5 6 2; 6 6 1 3], UnitHexagon), 6),
             (BindingRules([1 1 2 3; 2 1 3 3; 3 2 5 4; 5 2 8 4; 7 1 8 3; 6 1 7 3; 4 2 6 4; 1 2 4 4], UnitSquare), 8))

    for (rules, L) in rings
        ring = target(rules, L)
        @test nparticles(ring) == L
        @test nbonds(ring) == L # closed, so bonds == particles

        # Breaking one bond of a ring leaves a chain, which is still connected
        @test ncuts(ring; maxbonds=1) == 0
        # Every pair of bonds splits it into two arcs, and nothing larger is minimal
        @test ncuts(ring; maxbonds=2) == binomial(L, 2)
        @test ncuts(ring; maxbonds=Inf) == binomial(L, 2)
    end

    # A fully addressable 3x3 block of squares
    gridrules = BindingRules([1 1 2 3; 2 1 3 3; 4 1 5 3; 5 1 6 3; 7 1 8 3; 8 1 9 3;
                              1 2 4 4; 2 2 5 4; 3 2 6 4; 4 2 7 4; 5 2 8 4; 6 2 9 4], UnitSquare)
    grid = target(gridrules, 9)
    @test nparticles(grid) == 9
    @test nbonds(grid) == 12

    # Every bond of the grid lies on a four-cycle, so no single bond disconnects it.
    @test ncuts(grid; maxbonds=1) == 0
    # Only the four corner cells have as few as two bonds.
    @test ncuts(grid; maxbonds=2) == 4

    # Three-bond cuts: the four edge-centre cells, the eight corner-edge 2x1 blocks, and the two
    # outer rows and two outer columns, on top of the four corners above.
    @test ncuts(grid; maxbonds=3) == 4 + (4 + 8 + 4)

    # Larger cuts are not counted by hand, just pinned so the search cannot drift.
    @test ncuts(grid; maxbonds=4) == 37
    @test ncuts(grid; maxbonds=Inf) == 53

    # check that rates are non-negative
    rules = BindingRules([1 1 2 3], UnitTriangle)
    asys = AssemblySystem(rules, EntropyModel(TreeApproximation()))

    for bad in (-1.5, NaN, Inf, -Inf)
        @test_throws ArgumentError ReactionNetwork(asys; fwdkernel=Returns(bad))
        @test_throws ArgumentError ReactionNetwork(asys; fwdkernel=Returns(1.0), bwdkernel=Returns(bad))
    end

    net = ReactionNetwork(asys)
    for bad in (-1.5, NaN, Inf)
        @test_throws ArgumentError rate!(net; fwdkernel=Returns(bad))
    end
    
    # zero is a legitimate rate: it just deactivates the reaction
    @test !any(rate!(net; fwdkernel=Returns(0.0)).active)
    @test all(rate!(net; fwdkernel=Returns(1.0)).active)
end;
