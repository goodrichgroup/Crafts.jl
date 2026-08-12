@testset "environmentrelations" begin
    using NautyGraphs: labels
    dimerrules = BindingRules([1 1 2 1], UnitSquare)
    chainrules = BindingRules([1 1 2 3; 2 1 3 3], UnitSquare)
    central = PolygonParticleSpecies(3; labels=[1, 1, 1])
    outer = PolygonParticleSpecies(3)
    trianglerules = BindingRules([1 1 2 1; 1 2 2 1; 1 3 2 1], [central, outer])
    divalentrules = BindingRules([1 1 2 1; 1 3 2 1], UnitSquare)
    # one species forming unbounded chains: opposite sites 1 and 3 bind each other
    chainlike = BindingRules([1 1 1 3], UnitSquare)

    # normalize by a positive scalar only, so ray directions and facet orientations survive
    normdir(r) = r .// maximum(abs, r)
    dirset(rows) = Set(normdir(collect(r)) for r in eachrow(rows))
    dirset(rows::AbstractVector{<:AbstractVector}) = Set(normdir(r) for r in rows)

    # the 1-2 dimer system: four environments, one relation, and the exact outer cone
    rel = environmentrelations(dimerrules; depth=1)
    @test length(rel.envs) == 4
    @test size(rel.relations) == (1, 4)
    @test length(rel.bondclasses) == 1

    O = outercone(rel)
    @test dirset(O.rays) == dirset([[1, 0, 0], [0, 1, 0], [1, 1, 1]])
    @test dirset(O.facets) == dirset([[-1, 0, 1], [0, -1, 1], [0, 0, -1]])   # b <= n1, b <= n2, b >= 0

    # facet certificate: n1 - b reduces to the species-1 monomer coordinate under the relation
    c = [1, 0, -1]
    μ1 = environmentcounts(rel, Polyform(dimerrules, 1))
    @test any(s -> vec(c' * rel.projection) - s .* rel.relations[1, :] == μ1, (-1//2, 1//2))

    # every relation balances exactly on every explicit cluster, class by class
    for (rules, maxsize) in ((dimerrules, 4), (chainrules, 5), (trianglerules, 6),
                             (divalentrules, 5), (chainlike, 6))
        r1 = environmentrelations(rules; depth=1)
        for s in polygen(rules; maxsize)
            μ = environmentcounts(r1, s)
            @test sum(μ) == nparticles(s)
            @test all(iszero, r1.relations * μ)
        end
    end

    # the depth-2 union crop balances on asymmetric chain environments (the classic crop bug shows
    # up exactly here: cropping B_{k-1} of one endpoint instead of the union breaks the pairing)
    rel2 = environmentrelations(chainlike; depth=2)
    @test size(rel2.relations, 1) >= 2
    for s in polygen(chainlike; maxsize=6)
        μ = environmentcounts(rel2, s)
        @test all(iszero, rel2.relations * μ)
    end

    # the projection reproduces the composition, with bond counts grouped by symmetry class
    sitelabel(rules, loc) = labels(graphrep(Roly.species(rules, loc[1])))[first(
        bindingsites(Roly.species(rules, loc[1]), loc[2]).vertices)]
    function symcomposition(rel, rules, s)
        comp = composition(s)
        ns = nspecies(rules)
        sym = zeros(Rational{Int}, ns + length(rel.bondclasses))
        sym[1:ns] .= comp[1:ns]
        for (β, (locs1, locs2)) in enumerate(Roly.bonded_sites(rules))
            l1, l2 = first(locs1), first(locs2)
            cls = minmax((l1[1], sitelabel(rules, l1)), (l2[1], sitelabel(rules, l2)))
            sym[ns + findfirst(==(cls), rel.bondclasses)] += comp[ns + β]
        end
        return sym
    end
    for (rules, maxsize) in ((chainrules, 5), (trianglerules, 6), (divalentrules, 5))
        r1 = environmentrelations(rules; depth=1)
        for s in polygen(rules; maxsize)
            @test r1.projection * environmentcounts(r1, s) == symcomposition(r1, rules, s)
        end
    end
    # the triangle system's three bond types collapse into one symmetry class
    @test length(environmentrelations(trianglerules; depth=1).bondclasses) == 1

    # the two bond-counting conventions agree on every feasible count vector: half at both ends
    # versus a full count at the designated endpoint of each relation pair
    r1 = environmentrelations(chainrules; depth=1)
    designated = Dict(p[1] => 1//1 for p in r1.rowpairs)
    altbond(w) = sum((e == reverse(e) ? 1//2 : get(designated, e, 0//1)
                      for e in bondenvironments(w)); init=0//1)
    alt = [altbond(w) for w in r1.envs]
    for s in polygen(chainrules; maxsize=5)
        μ = environmentcounts(r1, s)
        @test alt' * μ == sum(r1.projection[(nspecies(chainrules) + 1):end, :] * μ)
    end

    # both remaining fixtures close at depth 1: the outer cone equals the hull of the cluster
    # compositions, with the fully assembled structure as an extreme ray
    Otri = outercone(environmentrelations(trianglerules; depth=1))
    @test dirset(Otri.rays) == dirset([[1, 0, 0], [0, 1, 0], [1, 3, 3]])

    rdiv = environmentrelations(divalentrules; depth=1)
    Odiv = outercone(rdiv)
    divcomps = Set(normdir(rdiv.projection * environmentcounts(rdiv, s))
                   for s in polygen(divalentrules; maxsize=5))
    @test dirset(Odiv.rays) == divcomps
    @test normdir([1, 2, 1, 1]) in dirset(Odiv.rays)   # the 2-1-2 trimer is forced to be extreme

    # exact LP solves: bounded, unbounded (with an improving ray), infeasible
    st, x, _ = solvelp([1 0; 0 1; -1 0; 0 -1], [1, 2, 0, 0], [1, 1])
    @test st == :optimal && x == [1, 2]
    st, _, rays = solvelp([-1 0; 0 -1], [0, 0], [1, 1])
    @test st == :unbounded && any(r -> sum(r) > 0, eachrow(rays))
    st, _, _ = solvelp(reshape([1, -1], 2, 1), [-1, 0], [1])
    @test st == :infeasible

    # gift-wrapping the octant through P recovers the projected quadrant
    P = [1 1 0; 0 0 1]
    Fq, Rq = projectcone(P, [-1 0 0; 0 -1 0; 0 0 -1])
    @test dirset(Rq) == dirset([[1, 0], [0, 1]])
    @test dirset(Fq) == dirset([[-1, 0], [0, -1]])

    # the LP route and the ray-enumeration route produce the same cone on every fixture
    for (rules, depth) in ((dimerrules, 1), (chainlike, 2), (divalentrules, 2), (trianglerules, 1))
        rel = environmentrelations(rules; depth)
        Or = outercone(rel; method=:rays)
        Op = outercone(rel; method=:project)
        @test dirset(Op.facets) == dirset(Or.facets)
        @test dirset(Op.rays) == dirset(Or.rays)
    end

    # the square-lattice system: at depth 2 the feasible cone 𝒦 has too many rays to enumerate,
    # but the projection is just the box 0 <= b1, b2 <= n — chains, lattice and monomer as rays
    squarerules = BindingRules([1 2 1 4; 1 3 1 1], UnitSquare)
    boxrays = dirset([[1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1]])
    Osq1 = outercone(environmentrelations(squarerules; depth=1))
    @test dirset(Osq1.rays) == boxrays
    relsq2 = environmentrelations(squarerules; depth=2)
    Osq2 = outercone(relsq2)   # :auto must pick :project here
    @test dirset(Osq2.rays) == boxrays
    @test dirset(Osq2.facets) == dirset(Osq1.facets)
    for s in polygen(squarerules; maxsize=5)
        m = relsq2.projection * environmentcounts(relsq2, s)
        @test all(x -> x <= 0, Osq2.facets * m)
    end

    # the float JuMP oracle reproduces the exact cones, with rationalization snapping exactly
    for (rules, depth) in ((chainrules, 1), (divalentrules, 2))
        rel = environmentrelations(rules; depth)
        Oe = outercone(rel; method=:project)
        Of = outercone(rel; optimizer=HiGHS.Optimizer)
        @test dirset(Of.facets) == dirset(Oe.facets)
        @test dirset(Of.rays) == dirset(Oe.rays)
    end
    @test dirset(outercone(relsq2; optimizer=HiGHS.Optimizer).rays) == boxrays

    # refinement strategy: the fiber support restricts the depth-2 test losslessly, so the
    # restricted survival verdict must match the full depth-2 relations on every depth-1 ray
    relsq1 = environmentrelations(squarerules; depth=1)
    for m in ([1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1])
        U = fibersupport(relsq1, m)
        @test !isempty(U)
        refined = refinementrelations(relsq1, U)
        @test refined.bondclasses == relsq1.bondclasses
        keptset = Set(relsq1.envs[U])
        @test all(crop(e, 1) in keptset for e in refined.envs)
        @test length(refined.envs) <= length(relsq2.envs)
        @test (realizablecounts(refined, m) !== nothing) ==
              (realizablecounts(relsq2, m) !== nothing)
    end
    # infeasible directions report cleanly
    @test realizablecounts(relsq1, [0, 1, 0]) === nothing
    @test isempty(fibersupport(relsq1, [0, 1, 0]))

    # a feasible count vector from the fiber balances and projects correctly
    μf = realizablecounts(relsq1, [1, 1, 1])
    @test μf !== nothing && all(>=(0), μf)
    @test all(iszero, relsq1.relations * μf) && relsq1.projection * μf == [1, 1, 1]

    # finite witnesses: every divalent ray is realized by an actual structure whose symmetrized
    # composition is proportional to it; the chainlike limit ray has no finite witness
    for r in eachrow(Odiv.rays)
        m = normdir(collect(r))
        w = findraywitness(rdiv, m; maxsize=6)
        @test w !== nothing
        @test normdir(rdiv.projection * environmentcounts(rdiv, w)) == m
    end
    relchain = environmentrelations(chainlike; depth=1)
    @test findraywitness(relchain, [1, 1]; maxscale=3, maxsize=10) === nothing

    # periodic witnesses certify the limit rays that no finite structure can: the witness is a
    # unit cell whose per-cell composition (bulk bonds plus the bonds the tiling closes) lands
    # on `m`. The chain ray is realized by the bare monomer; the square-lattice system realizes
    # its two chain rays as partial tilings and the full lattice ray as the complete one.
    wchain = findraywitness(relchain, [1, 1]; maxsize=4, periodic=true)
    @test wchain !== nothing && nparticles(wchain) == 1
    for m in ([1, 1, 0], [1, 0, 1], [1, 1, 1])
        @test findraywitness(relsq1, m; maxsize=4, periodic=true) !== nothing
    end
    # periodic search does not spuriously witness infeasible directions
    @test findraywitness(relsq1, [1, 2, 0]; maxsize=3, periodic=true) === nothing

    # the full certification loop: every fixture ray is witnessed, so C = O₁ is certified both
    # ways; the square lattice needs periodic witnesses for its three limit rays
    cert = certifyrays(rdiv)
    @test cert.certified && all(v -> v.status === :finite, cert.verdicts)
    certsq = certifyrays(relsq1)
    @test certsq.certified
    @test count(v -> v.status === :periodic, certsq.verdicts) == 3
    @test all(v -> v.witness !== nothing, certsq.verdicts)
    # without the periodic mode the chain ray survives refinement and stays undetermined
    certnp = certifyrays(relchain; periodic=false)
    @test !certnp.certified
    @test any(v -> v.status === :undetermined, certnp.verdicts)
    @test any(v -> v.status === :finite, certnp.verdicts)

    # seeding the gift-wrap with certified rays changes nothing about the resulting cone
    Wd = reduce(vcat, (permutedims(v.ray) for v in cert.verdicts))
    Oseed = outercone(rdiv; method=:project, seeds=Wd)
    @test dirset(Oseed.rays) == dirset(Odiv.rays)
    @test dirset(Oseed.facets) == dirset(Odiv.facets)
end
