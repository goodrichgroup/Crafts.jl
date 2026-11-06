const ExactNum = Union{Integer,Rational}

# The exact polyhedral primitives below are implemented by the Polymake package extension; the
# backend value only routes dispatch there.
struct NoBackend end
struct PolymakeBackend end

_backend() = Base.get_extension(@__MODULE__, :CraftsPolymakeExt) === nothing ? NoBackend() : PolymakeBackend()

_nopolymake() = throw(ArgumentError("exact polyhedral computation requires Polymake; run `using Polymake` first"))
_extremerays(::NoBackend, args...) = _nopolymake()
_facetsof(::NoBackend, args...) = _nopolymake()
_removeredundancy(::NoBackend, args...) = _nopolymake()
_solvelp(::NoBackend, args...) = _nopolymake()
_facelattice(::NoBackend, args...) = _nopolymake()

function _checkrows(A, b)
    size(A, 1) == length(b) ||
        throw(DimensionMismatch("`A` has $(size(A, 1)) rows but `b` has $(length(b))"))
    return nothing
end

"""
    extremerays(A, b; linearity=Int[])

Convert the H-representation `A*x <= b` into its extreme rays and vertices.

  - `linearity`: indices of rows of `A` that hold with equality
  - returns `(rays, vertices)`, both matrices with one generator per row

Requires Polymake to be loaded.
"""
function extremerays(A::AbstractMatrix{<:ExactNum}, b::AbstractVector{<:ExactNum}; linearity=Int[])
    _checkrows(A, b)
    return _extremerays(_backend(), A, b, linearity)
end

"""
    foreachray(f, A, b; linearity=Int[])

Apply `f(generator, isray)` to each generator of `A*x <= b`.

The generators are computed by [`extremerays`](@ref) first; this only iterates over them. Each
`generator` is a fresh vector and is safe to keep.

  - `linearity`: indices of rows of `A` that hold with equality
"""
function foreachray(f, A::AbstractMatrix{<:ExactNum}, b::AbstractVector{<:ExactNum}; linearity=Int[])
    rays, verts = extremerays(A, b; linearity)
    for r in eachrow(rays)
        f(collect(r), true)
    end
    for v in eachrow(verts)
        f(collect(v), false)
    end
    return nothing
end

"""
    facetsof(V; rays=true, incidence=true)

Convert the V-representation `V` into facets, with the incidence of generators on each facet.

  - `V`: one generator per row
  - `rays`: whether the rows are rays of a cone rather than vertices of a polytope; may also be
    one flag per row for a mixed representation
  - `incidence`: whether to compute the incidence at all; `false` returns an empty matrix in its
    place, which is what a caller that only wants the facets should pass
  - returns `(A, b, incidence)` with `A*x <= b`, and `incidence[i, j]` true when generator `j` lies
    on facet `i`

Only the facets are returned: the equations of a lower-dimensional cone are dropped. Incidence
costs one exact rational dot product per (facet, generator) pair, which dominates the conversion
itself on the large transient hulls of a gift-wrap; hence the switch. Requires Polymake to be
loaded.
"""
function facetsof(V::AbstractMatrix{<:ExactNum}; rays=true, incidence::Bool=true)
    m = size(V, 1)
    flags = rays isa Bool ? fill(rays, m) : collect(Bool, rays)
    length(flags) == m ||
        throw(DimensionMismatch("`rays` has $(length(flags)) entries but `V` has $m rows"))
    return _facetsof(_backend(), V, flags, incidence)
end

# Both representations are exact rationals, so testing `A*v == b` directly is cheaper to trust
# than any backend's own incidence report, whose generator order need not match the input.
function _incidence(A, b, V)
    inc = falses(size(A, 1), size(V, 1))
    for i in axes(A, 1), j in axes(V, 1)
        inc[i, j] = sum(A[i, k] * V[j, k] for k in axes(V, 2)) == b[i]
    end
    return inc
end

"""
    removeredundancy(A, b; linearity=Int[])

Drop the redundant inequalities from the H-representation `A*x <= b`.

  - `linearity`: indices of rows of `A` that hold with equality; these are always kept
  - returns `(A, b, kept)`, where `kept` indexes the surviving rows of the input

Of several equivalent rows exactly one survives, the first. Requires Polymake to be loaded.
"""
function removeredundancy(A::AbstractMatrix{<:ExactNum}, b::AbstractVector{<:ExactNum}; linearity=Int[])
    _checkrows(A, b)
    kept = _removeredundancy(_backend(), A, b, linearity)
    return A[kept, :], b[kept], kept
end

"""
    solvelp(A, b, w; linearity=Int[], maximize=true)

Solve `max wᵀx` (or min) subject to `A*x <= b` exactly.

  - `linearity`: indices of rows of `A` that hold with equality
  - returns `(status, x, rays)`: `status` is `:optimal`, `:unbounded` or `:infeasible`; `x` is an
    optimal vertex, or `nothing` otherwise; on `:unbounded`, `rays` holds one improving ray of
    the feasible region as its only row

The feasible region must be pointed. Requires Polymake to be loaded.
"""
function solvelp(A::AbstractMatrix{<:ExactNum}, b::AbstractVector{<:ExactNum},
                 w::AbstractVector{<:ExactNum}; linearity=Int[], maximize=true)
    _checkrows(A, b)
    length(w) == size(A, 2) ||
        throw(DimensionMismatch("`w` has $(length(w)) entries but `A` has $(size(A, 2)) columns"))
    return _solvelp(_backend(), A, b, w, linearity, maximize)
end

_tomatrix(rows, ncols) = isempty(rows) ? zeros(Rational{BigInt}, 0, ncols) :
                         permutedims(reduce(hcat, rows))

# Gauss-Jordan over the rationals, in place; returns the pivot columns, so the rank.
function _pivotcols!(M::AbstractMatrix{<:Rational})
    piv = Int[]
    r = 1
    for c in axes(M, 2)
        r > size(M, 1) && break
        i = findnext(!iszero, view(M, :, c), r)
        i === nothing && continue
        if i != r
            M[r, :], M[i, :] = M[i, :], M[r, :]
        end
        M[r, :] ./= M[r, c]
        for j in axes(M, 1)
            j == r || iszero(M[j, c]) || (M[j, :] .-= M[j, c] .* M[r, :])
        end
        push!(piv, c)
        r += 1
    end
    return piv
end

# A direction orthogonal to the row space of `G`, or `nothing` when `G` has full column rank.
function _orthdir(G::AbstractMatrix{<:Rational})
    M = copy(G)
    piv = _pivotcols!(M)
    f = findfirst(!in(piv), 1:size(G, 2))
    f === nothing && return nothing
    x = zeros(Rational{BigInt}, size(G, 2))
    x[f] = 1
    for (k, c) in enumerate(piv)
        x[c] = -M[k, f]
    end
    return x
end

# Scale by a positive factor to a canonical representative of the direction, for dedup keys.
_normdir(v) = v ./ maximum(abs, v)

# The exact support-LP oracle: `oracle(c)` returns `(valid, points)` where `valid` means
# `sup ⟨c, P*x⟩ = 0` over the source cone, and `points` are projected improving generators.
function _lporacle(::Nothing, P, A, linearity)
    b = zeros(Rational{Int}, size(A, 1))
    return function (c)
        status, _, rays = solvelp(A, b, transpose(P) * c; linearity)
        status === :optimal && return true, Vector{Rational{BigInt}}[]
        return false, [P * collect(r) for r in eachrow(rays)]
    end
end

_lporacle(optimizer, P, A, linearity) =
    throw(ArgumentError("passing an `optimizer` requires JuMP to be loaded, e.g. `using JuMP, HiGHS`"))

"""
    projectcone(P, A; linearity=Int[], seeds=nothing, optimizer=nothing)

Compute the facets and extreme rays of the projected cone `{P*x : A*x <= 0}` by LP gift-wrapping,
without ever enumerating the rays of the source cone.

  - `linearity`: indices of rows of `A` that hold with equality
  - `seeds`: optional known points of the projected cone, one per row, to start the hull from
  - `optimizer`: a JuMP-compatible optimizer (e.g. `HiGHS.Optimizer`; requires JuMP loaded) to
    solve the support LPs in floating point, with results rationalized; much faster than the
    default exact solves, but the face certificates are no longer exact
  - returns `(facets, rays)`, one inequality `f ⋅ y <= 0` and one ray per row

The projected cone must be full-dimensional. The cost is one LP per face certificate plus one per
generator batch, so it scales with the size of the projected cone, not the source cone. Requires
Polymake to be loaded.
"""
function projectcone(P::AbstractMatrix{<:ExactNum}, A::AbstractMatrix{<:ExactNum};
                     linearity=Int[], seeds=nothing, optimizer=nothing)
    d, n = size(P)
    size(A, 2) == n ||
        throw(DimensionMismatch("`P` has $n columns but `A` has $(size(A, 2))"))
    oracle = _lporacle(optimizer, P, A, linearity)

    gens = Vector{Rational{BigInt}}[]
    seen = Set{Vector{Rational{BigInt}}}()
    certified = Set{Vector{Rational{BigInt}}}()
    function addgen(v)
        iszero(v) && return false
        key = _normdir(v)
        key in seen && return false
        push!(seen, key)
        push!(gens, key)
        return true
    end
    seeds === nothing || foreach(s -> addgen(collect(Rational{BigInt}, s)), eachrow(seeds))

    # Either `c` is valid for the whole projection (`:valid`), or the oracle hands back improving
    # generators that grow the hull (`:grew`), or it disputes `c` without one (`:stuck`). A
    # dispute only counts as growth when a returned point strictly violates `c` in exact
    # arithmetic: a float oracle can dispute within tolerance yet return a point that does not
    # actually beat the hyperplane, and admitting such points would recompute the identical hull
    # and re-probe the identical facet forever.
    function probe(c)
        valid, points = oracle(c)
        valid && return :valid
        grew = false
        for y in points
            sum(c .* y) > 0 && (grew |= addgen(y))
        end
        return grew ? :grew : :stuck
    end

    # Probe each candidate facet not already settled in `settled`; recompute the hull as soon as
    # it grows. Probing the remaining facets of a stale hull is counterproductive, their junk
    # generators balloon the hull far beyond what the saved recomputation costs.
    function sweep!(settled, candidates, key)
        grew = false
        for a in eachrow(candidates)
            f = collect(a)
            k = key(f)
            k in settled && continue
            r = probe(f)
            if r === :grew
                grew = true
                break
            end
            r === :stuck &&
                @warn "the oracle disputed a facet without providing a new generator; accepting it" maxlog=3
            push!(settled, k)
        end
        return grew
    end

    # Bootstrap to full dimension: probe both sides of a direction orthogonal to the current hull.
    while (c = _orthdir(_tomatrix(gens, d))) !== nothing
        probe(c) === :grew && continue
        probe(-c) === :grew && continue
        throw(ArgumentError("the projected cone is not full-dimensional"))
    end

    # Grow the hull until every facet of it is certified valid for the projection. The candidate
    # facets always come from the exact hull: float hull candidates were tried and are unstable
    # against the support oracle. Over a cone, facet validity is discontinuous in the direction,
    # so noisy candidates read as violated and flood the hull with junk generators, however the
    # tolerance is arranged.
    while true
        F, _, _ = facetsof(_tomatrix(gens, d); rays=true, incidence=false)
        sweep!(certified, F, _normdir) && continue
        rays, _ = extremerays(F, zeros(Rational{Int}, size(F, 1)))
        return F, rays
    end
end
