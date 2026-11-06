module CraftsPolymakeExt

using Crafts, Polymake
using Crafts: PolymakeBackend

const pm = Polymake

# polymake solves rational LPs with cdd by default, which is orders of magnitude slower than
# SoPlex on the support LPs of the gift-wrap (a 1185×1057 fiber system: cdd 555 s, lrs 46 s,
# polymake's own simplex 2.8 s, SoPlex 1.2 s). The preference is global to the polymake session.
function __init__()
    app = pm.get_current_app()
    pm.shell_execute("application \"polytope\"; prefer \"soplex.simplex\";")
    pm.application(app)
    return nothing
end

_pmmatrix(M) = pm.Matrix{pm.Rational}(Matrix{Rational{BigInt}}(M))
_pmvector(v) = pm.Vector{pm.Rational}(Vector{Rational{BigInt}}(v))
_jlmatrix(M) = Matrix{Rational{BigInt}}(M)
_jlvector(v) = Vector{Rational{BigInt}}(v)

# polymake writes a half-space as `h ⋅ [1; x] >= 0`, so `A*x <= b` becomes the rows `[b -A]`.
# Returns the inequality rows and the equation rows.
function _hrep(A, b, linearity)
    islin = falses(size(A, 1))
    islin[linearity] .= true
    H = hcat(Rational{BigInt}.(b), -Rational{BigInt}.(A))
    return H[.!islin, :], H[islin, :]
end

_nonzerorows(H) = [any(!iszero, @view H[i, 2:end]) for i in axes(H, 1)]

# polymake refuses all-zero rows. A zero inequality `0 >= -b` or equation `0 = b` either holds
# trivially or makes the system infeasible, so such rows are dropped, or the whole system is
# replaced by one unsatisfiable inequality.
function _cleanrows(ineq, eqs)
    d = size(ineq, 2)
    ki = _nonzerorows(ineq)
    ke = _nonzerorows(eqs)
    contradiction = any(i -> !ki[i] && ineq[i, 1] < 0, axes(ineq, 1)) ||
                    any(i -> !ke[i] && !iszero(eqs[i, 1]), axes(eqs, 1))
    contradiction && return reshape([-one(Rational{BigInt}); zeros(Rational{BigInt}, d - 1)], 1, d),
                            zeros(Rational{BigInt}, 0, d)
    return ineq[ki, :], eqs[ke, :]
end

function _polytope(ineq, eqs)
    ineq, eqs = _cleanrows(ineq, eqs)
    return pm.polytope.Polytope{pm.Rational}(INEQUALITIES=_pmmatrix(ineq), EQUATIONS=_pmmatrix(eqs))
end

function _linearprogram(P, obj)
    P.LP = pm.polytope.LinearProgram{pm.Rational}(LINEAR_OBJECTIVE=_pmvector(obj))
    return P.LP
end

function Crafts._extremerays(::PolymakeBackend, A, b, linearity)
    d = size(A, 2)
    ineq, eqs = _hrep(A, b, linearity)
    if all(iszero, b)
        # a homogeneous zero row holds trivially, and polymake refuses it
        C = pm.polytope.Cone{pm.Rational}(INEQUALITIES=_pmmatrix(ineq[_nonzerorows(ineq), 2:end]),
                                          EQUATIONS=_pmmatrix(eqs[_nonzerorows(eqs), 2:end]))
        return _jlmatrix(C.RAYS), zeros(Rational{BigInt}, 0, d)
    end
    P = _polytope(ineq, eqs)
    P.FEASIBLE || return zeros(Rational{BigInt}, 0, d), zeros(Rational{BigInt}, 0, d)
    V = _jlmatrix(P.VERTICES)
    isray = iszero.(V[:, 1])
    rays = V[isray, 2:end]
    verts = V[.!isray, 2:end] ./ V[.!isray, 1]
    return rays, verts
end

function Crafts._facetsof(::PolymakeBackend, V, flags, incidence)
    d = size(V, 2)
    if all(flags)
        C = pm.polytope.Cone{pm.Rational}(INPUT_RAYS=_pmmatrix(V))
        A = -_jlmatrix(C.FACETS)
        b = zeros(Rational{BigInt}, size(A, 1))
    else
        # vertices get the homogenizing coordinate 1, rays 0; the far face `x₀ >= 0` of an
        # unbounded polyhedron comes back as a facet and is dropped
        pts = hcat(Rational{BigInt}.(.!flags), Rational{BigInt}.(V))
        P = pm.polytope.Polytope{pm.Rational}(POINTS=_pmmatrix(pts))
        F = _jlmatrix(P.FACETS)
        proper = [any(!iszero, @view F[i, 2:end]) for i in axes(F, 1)]
        A = -F[proper, 2:end]
        b = F[proper, 1]
    end
    isempty(A) && (A = zeros(Rational{BigInt}, 0, d))
    return A, b, incidence ? Crafts._incidence(A, b, V) : falses(size(A, 1), 0)
end

function Crafts._removeredundancy(::PolymakeBackend, A, b, linearity)
    m = size(A, 1)
    islin = falses(m)
    islin[linearity] .= true
    H = hcat(Rational{BigInt}.(b), -Rational{BigInt}.(A))
    keep = trues(m)

    P = _polytope(H[.!islin, :], H[islin, :])
    if P.FULL_DIM
        # Facet normals of a full-dimensional polyhedron are unique up to positive scaling, so a
        # row is irredundant iff it is a facet, and the first row on each facet is the survivor.
        F = _jlmatrix(P.FACETS)
        facets = Set(Crafts._normdir(f) for f in eachrow(F) if any(!iszero, @view f[2:end]))
        for i in 1:m
            islin[i] && continue
            if any(!iszero, @view H[i, 2:end])
                key = Crafts._normdir(H[i, :])
                keep[i] = key in facets && pop!(facets, key) !== nothing
            else
                keep[i] = false
            end
        end
    else
        # Within an affine subspace a row is redundant iff maximizing it over the other rows
        # cannot exceed its bound; one LP per row, against the rows still kept.
        for i in 1:m
            islin[i] && continue
            others = [j for j in 1:m if keep[j] && !islin[j] && j != i]
            bound = permutedims(vcat(b[i] + 1, -Rational{BigInt}.(A[i, :])))
            Pi = _polytope(vcat(H[others, :], bound), H[islin, :])
            lp = _linearprogram(Pi, vcat(0, A[i, :]))
            v = lp.MAXIMAL_VALUE
            keep[i] = v === nothing || Rational{BigInt}(v) > b[i]
        end
    end
    return findall(keep)
end

function Crafts._solvelp(::PolymakeBackend, A, b, w, linearity, maximize)
    d = size(A, 2)
    norays = zeros(Rational{BigInt}, 0, d)
    ineq, eqs = _hrep(A, b, linearity)
    P = _polytope(ineq, eqs)
    obj = maximize ? Rational{BigInt}.(w) : -Rational{BigInt}.(w)
    lp = _linearprogram(P, vcat(0, obj))
    value = lp.MAXIMAL_VALUE   # `nothing` on an infeasible system
    value === nothing && return :infeasible, nothing, norays
    if !isinf(value)
        x = lp.MAXIMAL_VERTEX
        x === nothing && throw(ArgumentError("the feasible region has lineality; `solvelp` needs a pointed one"))
        return :optimal, _jlvector(x)[2:end], norays
    end
    # An improving ray is a vertex of the recession cone sliced at `obj ⋅ r = 1`.
    rec = vcat(hcat(zeros(Rational{BigInt}, size(ineq, 1)), ineq[:, 2:end]),
               permutedims(vcat(one(Rational{BigInt}), -obj)))
    R = _polytope(rec, hcat(zeros(Rational{BigInt}, size(eqs, 1)), eqs[:, 2:end]))
    ray = _linearprogram(R, vcat(0, obj)).MAXIMAL_VERTEX
    ray === nothing && throw(ArgumentError("the feasible region has lineality; `solvelp` needs a pointed one"))
    return :unbounded, nothing, permutedims(_jlvector(ray)[2:end])
end

# Ranks of the Hasse diagram's nodes, by breadth-first search from its bottom; the diagram is
# graded, so every edge raises the rank by one. Node indices are polymake's, zero-based.
function _ranks(HD)
    n = HD.N_NODES
    G = HD.ADJACENCY
    rank = fill(-1, n)
    queue = [Int(HD.BOTTOM_NODE)]
    rank[queue[1] + 1] = 0
    while !isempty(queue)
        u = popfirst!(queue)
        for v in pm._outneighbors(G, u)
            rank[v + 1] == -1 || continue
            rank[v + 1] = rank[u + 1] + 1
            push!(queue, v)
        end
    end
    return rank
end

function Crafts._facelattice(::PolymakeBackend, inc::Matrix{Bool}, maxdim, covers, facetsets)
    nfacets, ngens = size(inc)
    # Generators on every facet belong to every face and facets through every generator bound
    # nothing; polymake's lattice is built for incidences without either, so they are set aside
    # and the shared generators are put back into each face afterwards.
    shared = vec(all(inc; dims=1))
    proper = .!vec(all(inc; dims=2))
    keptgens = findall(.!shared)
    IM = pm.IncidenceMatrix(inc[proper, keptgens])
    dim = Int(pm.call_function(:polytope, :dim_from_incidence, IM))
    # `lower_hasse_diagram` closes the truncated lattice with an artificial top, marked by the
    # generator index -1
    HD = maxdim > dim ? pm.call_function(:polytope, :hasse_diagram, IM, dim) :
         pm.call_function(:polytope, :lower_hasse_diagram, IM, maxdim)
    n = HD.N_NODES
    rank = _ranks(HD)

    facesets = Vector{BitVector}(undef, n)
    genuine = trues(n)
    for i in 1:n
        members = collect(HD.FACES[i])
        if any(<(0), members)
            genuine[i] = false
            continue
        end
        f = copy(shared)
        for g in members
            f[keptgens[g + 1]] = true
        end
        facesets[i] = f
    end

    nlevels = maximum(rank[genuine]) + 1
    levels = [BitVector[] for _ in 1:nlevels]
    duals = facetsets ? [BitVector[] for _ in 1:nlevels] : nothing
    position = zeros(Int, n)
    for i in 1:n
        genuine[i] || continue
        level = levels[rank[i] + 1]
        push!(level, facesets[i])
        position[i] = length(level)
        facetsets || continue
        # the facets containing a face are those containing every one of its generators
        push!(duals[rank[i] + 1], BitVector(all(j -> !facesets[i][j] || inc[r, j], 1:ngens) for r in 1:nfacets))
    end

    edges = nothing
    if covers
        edges = [Tuple{Int,Int}[] for _ in 1:(nlevels - 1)]
        G = HD.ADJACENCY
        for u in 0:(n - 1)
            genuine[u + 1] || continue
            for v in pm._outneighbors(G, u)
                genuine[v + 1] || continue
                push!(edges[rank[u + 1] + 1], (position[u + 1], position[v + 1]))
            end
        end
    end
    return levels, duals, edges
end

end
