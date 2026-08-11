# Floating-point support-hyperplane candidates via qhull, for the gift-wrap loop in projectcone.
# Candidate facets need no exactness — a wrong candidate is either certified by the LP oracle or
# grows the generator set — so the many transient hulls along the way can be computed in
# milliseconds here instead of seconds in exact lrs arithmetic.

# Facet normals of `cone(V)` (one ray generator per row) in floating point, one normal per row,
# oriented as `n ⋅ x <= 0` inside. The hull is taken over the origin plus the generators; facets
# through the origin are the cone facets, the compact hull's cap facets are dropped by their
# offset. Retries with joggled input when qhull hits a precision error.
function _qhullconefacets(V::AbstractMatrix{Float64}; offsettol=1e-8, joggle=false)
    d = size(V, 2)
    io = IOBuffer()
    println(io, d)
    println(io, size(V, 1) + 1)
    println(io, join(zeros(Int, d), ' '))
    for r in eachrow(V)
        println(io, join(r, ' '))
    end
    opts = joggle ? ["n", "QJ"] : ["n"]

    out = IOBuffer()
    quiet = pipeline(`$(qconvex()) $opts`; stderr=devnull)
    try
        run(pipeline(pipeline(IOBuffer(String(take!(io))), quiet); stdout=out))
    catch
        joggle && rethrow()
        return _qhullconefacets(V; offsettol, joggle=true)
    end

    normals = Vector{Float64}[]
    for (i, line) in enumerate(eachline(IOBuffer(String(take!(out)))))
        i <= 2 && continue   # header: dimension+1, facet count
        row = [parse(Float64, t) for t in eachsplit(line)]
        length(row) == d + 1 || continue
        abs(row[end]) <= offsettol && push!(normals, row[1:d])
    end
    return permutedims(reduce(hcat, normals))
end
