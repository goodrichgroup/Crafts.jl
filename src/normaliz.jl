# Exact H<->V cone conversions through Normaliz. Chosen over lrs for these because its
# pyramid-decomposition algorithms handle the degenerate hulls that stall reverse search, and it
# parallelizes. Normaliz only speaks files: each call writes its input into a fresh temporary
# directory, runs the binary, and parses the output files back. Rational entries are accepted
# verbatim; results come back as primitive integer vectors.
#
# Conventions: inequality rows mean `row ⋅ x >= 0`; in inhomogeneous data the constant is the
# LAST column, and `.ext` rows carry the dehomogenizing coordinate last (0 for a ray).

const _NMZ_EMPTY = zeros(Rational{BigInt}, 0, 0)

function _nmzmatrix(tokens, i=1)
    i > length(tokens) && return _NMZ_EMPTY, i
    m = parse(Int, tokens[i])
    n = parse(Int, tokens[i + 1])
    M = Matrix{Rational{BigInt}}(undef, m, n)
    k = i + 2
    for r in 1:m, c in 1:n
        M[r, c] = _parserational(tokens[k])
        k += 1
    end
    return M, k
end

# `.cst` holds a sequence of blocks, each a matrix followed by its name.
function _nmzblocks(tokens)
    blocks = Dict{String,Matrix{Rational{BigInt}}}()
    i = 1
    while i <= length(tokens)
        M, i = _nmzmatrix(tokens, i)
        blocks[String(tokens[i])] = M
        i += 1
    end
    return blocks
end

_nmztokens(path) = isfile(path) ? split(read(path, String)) : SubString{String}[]

function _runnmz(ambdim, blocks, goals)
    mktempdir() do dir
        base = joinpath(dir, "nmz")
        open(base * ".in", "w") do io
            println(io, "amb_space ", ambdim)
            for (name, M) in blocks
                size(M, 1) == 0 && continue
                println(io, name, ' ', size(M, 1))
                _writerows(io, eachrow(M))
            end
            foreach(g -> println(io, g), goals)
        end
        # a spawn is occasionally killed silently under load; the input is fine, so retry
        ok = false
        local msg
        for attempt in 1:3
            err = IOBuffer()
            ok = success(pipeline(`$(normaliz()) -a $base`; stdout=devnull, stderr=err))
            msg = strip(String(take!(err)))
            ok && break
        end
        if !ok
            kept = tempname() * "-normaliz-failed.in"
            cp(base * ".in", kept)
            error("normaliz failed 3 times", isempty(msg) ? " with no message (killed?)" : ": $msg",
                  "; input saved to $kept")
        end
        return _nmzmatrix(_nmztokens(base * ".ext"))[1], _nmzblocks(_nmztokens(base * ".cst"))
    end
end

"""
    extremerays(A, b; linearity=Int[])

Convert the H-representation `A*x <= b` into its extreme rays and vertices.

  - `linearity`: indices of rows of `A` that hold with equality
  - returns `(rays, vertices)`, both matrices with one generator per row
"""
function extremerays(A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum}; linearity=Int[])
    m, d = size(A)
    length(b) == m ||
        throw(DimensionMismatch("`A` has $m rows but `b` has $(length(b))"))
    islin = falses(m)
    islin[linearity] .= true

    if all(iszero, b)
        ext, _ = _runnmz(d, ["inequalities" => -A[.!islin, :], "equations" => A[islin, :]],
                         ["ExtremeRays"])
        return ext, zeros(Rational{BigInt}, 0, d)
    end
    ext, _ = _runnmz(d,
                     ["inhom_inequalities" => hcat(-A[.!islin, :], b[.!islin]),
                      "inhom_equations" => hcat(A[islin, :], -b[islin])],
                     ["VerticesOfPolyhedron", "ExtremeRays"])
    rays = [collect(r[1:d]) for r in eachrow(ext) if iszero(r[end])]
    verts = [collect(r[1:d] ./ r[end]) for r in eachrow(ext) if !iszero(r[end])]
    return _tomatrix(rays, d), _tomatrix(verts, d)
end

"""
    facetsof(V; rays=true)

Convert the V-representation `V` into facets, with the incidence of generators on each facet.

  - `V`: one generator per row
  - `rays`: whether the rows are rays of a cone rather than vertices of a polytope; may also be
    one flag per row for a mixed representation
  - returns `(A, b, incidence)` with `A*x <= b`, and `incidence[i, j]` true when generator `j` lies
    on facet `i`
"""
function facetsof(V::AbstractMatrix{<:LrsNum}; rays=true)
    m, d = size(V)
    flags = rays isa Bool ? fill(rays, m) : collect(Bool, rays)
    length(flags) == m ||
        throw(DimensionMismatch("`rays` has $(length(flags)) entries but `V` has $m rows"))

    if all(flags)
        _, cst = _runnmz(d, ["cone" => V], ["SupportHyperplanes"])
        ineq = get(cst, "inequalities", _NMZ_EMPTY)
        A = isempty(ineq) ? zeros(Rational{BigInt}, 0, d) : -ineq
        b = zeros(Rational{BigInt}, size(A, 1))
    else
        nverts = count(!, flags)
        _, cst = _runnmz(d,
                         ["vertices" => hcat(V[.!flags, :], ones(Int, nverts)),
                          "cone" => V[flags, :]],
                         ["SupportHyperplanes"])
        ineq = get(cst, "inequalities", _NMZ_EMPTY)
        A = isempty(ineq) ? zeros(Rational{BigInt}, 0, d) : -ineq[:, 1:d]
        b = isempty(ineq) ? Rational{BigInt}[] : ineq[:, d + 1]
    end
    return A, b, _incidence(A, b, V)
end
