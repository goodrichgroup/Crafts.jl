@testset "lrs" begin
    # triangle {x >= 0, y >= 0, x + y <= 1}, with a redundant scaled copy of the third row
    A = Rational{BigInt}[-1 0; 0 -1; 1 1; 2 2]
    b = Rational{BigInt}[0, 0, 1, 2]

    @testset "extremerays" begin
        rays, verts = extremerays(A, b)
        @test isempty(rays)
        @test size(verts) == (3, 2)
        @test sort(collect(eachrow(verts))) == [[0, 0], [0, 1], [1, 0]]
    end

    @testset "removeredundancy" begin
        _, _, kept = removeredundancy(A, b)
        @test length(kept) == 3
        # rows 3 and 4 are equivalent, so exactly one of them survives
        @test count(in(kept), (3, 4)) == 1
        @test 1 in kept && 2 in kept

        _, _, allkept = removeredundancy(A[1:3, :], b[1:3])
        @test allkept == [1, 2, 3]
    end

    @testset "facetsof" begin
        V = Rational{BigInt}[0 0; 1 0; 1 1; 0 1]
        Af, bf, inc = facetsof(V; rays=false)
        @test size(Af, 1) == 4
        @test size(inc) == (4, 4)
        # every edge of a square holds exactly two vertices, and every vertex two edges
        @test all(==(2), sum(inc; dims=2))
        @test all(==(2), sum(inc; dims=1))
        # incidence agrees with evaluating the constraints
        for i in axes(Af, 1), j in axes(V, 1)
            @test inc[i, j] == (sum(Af[i, k] * V[j, k] for k in axes(V, 2)) == bf[i])
        end
    end

    @testset "round trip" begin
        V = Rational{BigInt}[0 0; 1 0; 1 1; 0 1]
        Af, bf, _ = facetsof(V; rays=false)
        _, verts = extremerays(Af, bf)
        @test sort(collect(eachrow(verts))) == sort(collect(eachrow(V)))
    end

    @testset "foreachray" begin
        seen = Vector{Rational{BigInt}}[]
        flags = Bool[]
        foreachray(A, b) do g, isray
            push!(seen, g)
            push!(flags, isray)
        end
        _, verts = extremerays(A, b)
        @test sort(seen) == sort(collect(eachrow(verts)))
        @test !any(flags)

        # the buffer handed to the callback is reused internally, so the copies must be distinct
        @test allunique(objectid.(seen))
    end

    @testset "cone rays" begin
        # {x >= 0, y >= 0, z >= 0, x + y - z = 0}
        Ac = Rational{BigInt}[-1 0 0; 0 -1 0; 0 0 -1; 1 1 -1]
        bc = Rational{BigInt}[0, 0, 0, 0]
        rays, _ = extremerays(Ac, bc; linearity=[4])
        normed = sort([r ./ maximum(r) for r in eachrow(rays) if any(!iszero, r)])
        @test normed == [[0, 1, 1], [1, 0, 1]]
    end
end
