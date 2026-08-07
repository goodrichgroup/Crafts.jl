@testset "facelattice" begin
    cubeverts(d) = Rational{BigInt}[((i >> (k - 1)) & 1) for i in 0:(2^d - 1), k in 1:d]
    simplexverts(d) = collect(Rational{BigInt}[k == 0 ? 0 : (j == k) for j in 1:d, k in 0:d]')

    _binom(n, k) = factorial(big(n)) ÷ (factorial(big(k)) * factorial(big(n - k)))
    expected_cube(d) = [1; [2^(d - k) * _binom(d, k) for k in 0:d]]
    expected_simplex(d) = [1; [_binom(d + 1, k + 1) for k in 0:d]]

    incidenceof(V) = last(facetsof(V; rays=false))

    @testset "cube f-vectors" begin
        for d in 2:4
            L = facelattice(incidenceof(cubeverts(d)))
            @test fvector(L) == expected_cube(d)
        end
    end

    @testset "simplex f-vectors" begin
        for d in 2:5
            L = facelattice(incidenceof(simplexverts(d)))
            @test fvector(L) == expected_simplex(d)
        end
    end

    inc4 = incidenceof(cubeverts(4))
    full = fvector(facelattice(inc4))

    @testset "truncation" begin
        for k in 0:4
            f = fvector(facelattice(inc4; maxdim=k))
            @test f == full[1:min(k + 1, length(full))]
        end
    end

    @testset "duality" begin
        @test fvector(facelattice(inc4; dual=true)) == reverse(full)
    end

    @testset "faces and structure" begin
        L = facelattice(incidenceof(cubeverts(3)); covers=true, facetsets=true)
        # a 3-cube: 8 vertices, 12 edges each on 2 vertices, 6 square facets each on 4
        @test length(faces(L, 0)) == 1 && !any(faces(L, 0)[1])
        @test all(f -> count(f) == 1, faces(L, 1))
        @test all(f -> count(f) == 2, faces(L, 2))
        @test all(f -> count(f) == 4, faces(L, 3))
        @test count(faces(L, 4)[1]) == 8

        # every edge covers two vertices, every square four edges
        @test length(covers(L, 1)) == 2 * 12
        @test length(covers(L, 2)) == 4 * 6

        # each vertex of a 3-cube lies on 3 facets, each edge on 2, each facet on 1
        @test all(f -> count(f) == 3, L.dualsets[2])
        @test all(f -> count(f) == 2, L.dualsets[3])
        @test all(f -> count(f) == 1, L.dualsets[4])
    end

    @testset "covers not retained" begin
        L = facelattice(incidenceof(simplexverts(3)))
        @test_throws ArgumentError covers(L, 1)
    end

    @testset "from an H-representation" begin
        # the unit square as inequalities
        A = Rational{BigInt}[-1 0; 0 -1; 1 0; 0 1]
        b = Rational{BigInt}[0, 0, 1, 1]
        @test fvector(facelattice(A, b)) == [1, 4, 4, 1]
    end
end
