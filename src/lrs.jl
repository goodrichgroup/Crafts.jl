const LrsNum = Union{Integer,Rational}

_printnum(io, x) = x isa Rational ? print(io, numerator(x), '/', denominator(x)) : print(io, x)

function _writerows(io, rows)
    for r in rows
        for x in r
            print(io, ' ')
            _printnum(io, x)
        end
        println(io)
    end
end

# lrs H-representation: each row is [b -A...], meaning b - A*x >= 0.
function _ine(A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum}; linearity=Int[],
              options=String[])
    size(A, 1) == length(b) ||
        throw(DimensionMismatch("`A` has $(size(A, 1)) rows but `b` has $(length(b))"))
    io = IOBuffer()
    println(io, "crafts\nH-representation")
    isempty(linearity) || println(io, "linearity ", length(linearity), " ", join(linearity, " "))
    println(io, "begin")
    println(io, size(A, 1), " ", size(A, 2) + 1, " rational")
    _writerows(io, (vcat(b[i], -A[i, :]) for i in axes(A, 1)))
    println(io, "end")
    for o in options
        println(io, o)
    end
    return String(take!(io))
end

function _parserational(s::AbstractString)
    i = findfirst('/', s)
    i === nothing && return Rational{BigInt}(parse(BigInt, s))
    return Rational{BigInt}(parse(BigInt, SubString(s, 1, i - 1)),
                            parse(BigInt, SubString(s, i + 1)))
end

# Cheaper than a regex and allocation-free: a data row is digits, minus signs, slashes and blanks.
function _isnumericrow(s::AbstractString)
    isempty(s) && return false
    for c in s
        (isdigit(c) || c == '-' || c == '/' || c == ' ' || c == '\t') || return false
    end
    return true
end

function _parserow!(buf::Vector{Rational{BigInt}}, s::AbstractString)
    empty!(buf)
    for t in eachsplit(s)
        push!(buf, _parserational(t))
    end
    return buf
end

# Runs lrs and hands back one row at a time, so neither the raw output nor the parsed rows are ever
# held in full. `onrow` receives a buffer that is reused between calls; copy it to keep it.
#
# lrs writes the answer as bare numeric lines terminated by `end`, and only afterwards reports things
# like which rows were redundant. Those trailing lines go to `ontrailer`, and must not be mistaken
# for data: the redundant-row indices look exactly like a data row.
function _streamlrs(input::AbstractString; onrow=nothing, ontrailer=nothing)
    buf = Rational{BigInt}[]
    # lrs writes a banner and progress notes to stderr; without this they land on the user's terminal.
    quiet = pipeline(`$(lrs())`; stderr=devnull)
    open(pipeline(IOBuffer(input), quiet)) do io
        trailing = false
        for line in eachline(io)
            s = strip(line)
            if trailing
                ontrailer === nothing || ontrailer(s)
                continue
            end
            (isempty(s) || startswith(s, '*') || startswith(s, "F#")) && continue
            if s == "end"
                trailing = true
                continue
            end
            _isnumericrow(s) || continue
            onrow === nothing || onrow(_parserow!(buf, s))
        end
    end
    return nothing
end

_tomatrix(rows, ncols) = isempty(rows) ? zeros(Rational{BigInt}, 0, ncols) :
                         permutedims(reduce(hcat, rows))

"""
    foreachray(f, A, b; linearity=Int[])

Apply `f(generator, isray)` to each generator of `A*x <= b` as lrs produces it.

Generators are never accumulated, so this needs far less memory than [`extremerays`](@ref) on systems
with many of them. `generator` is a fresh vector each call and is safe to keep.

  - `linearity`: indices of rows of `A` that hold with equality
"""
function foreachray(f, A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum}; linearity=Int[])
    _streamlrs(_ine(A, b; linearity); onrow=r -> f(r[2:end], iszero(r[1])))
    return nothing
end

# lrs can report incidences itself, but its `vertices/rays` lines mark membership with a `*` whose
# meaning shifts on degenerate input. Both representations are exact rationals here, so testing
# `A*v == b` directly is cheaper to trust.
function _incidence(A, b, V)
    inc = falses(size(A, 1), size(V, 1))
    for i in axes(A, 1), j in axes(V, 1)
        inc[i, j] = sum(A[i, k] * V[j, k] for k in axes(V, 2)) == b[i]
    end
    return inc
end

"""
    solvelp(A, b, w; linearity=Int[], maximize=true)

Solve `max wᵀx` (or min) subject to `A*x <= b` exactly.

  - `linearity`: indices of rows of `A` that hold with equality
  - returns `(status, x, rays)`: `status` is `:optimal`, `:unbounded` or `:infeasible`; `x` is the
    terminal vertex; on `:unbounded`, `rays` holds feasible rays at `x`, at least one improving
"""
function solvelp(A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum},
                 w::AbstractVector{<:LrsNum}; linearity=Int[], maximize=true)
    io = IOBuffer()
    print(io, maximize ? "maximize" : "minimize", " 0")
    for x in w
        print(io, ' ')
        _printnum(io, x)
    end
    x = nothing
    rays = Vector{Rational{BigInt}}[]
    _streamlrs(_ine(A, b; linearity, options=[String(take!(io)), "lponly"]);
               onrow=r -> iszero(r[1]) ? push!(rays, r[2:end]) : (x = r[2:end]))
    status = x === nothing ? :infeasible : isempty(rays) ? :optimal : :unbounded
    return status, x, _tomatrix(rays, size(A, 2))
end

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
    solve the support LPs in floating point, with results rationalized — much faster than the
    default exact lrs solves, but the face certificates are no longer exact
  - returns `(facets, rays)`, one inequality `f ⋅ y <= 0` and one ray per row

The projected cone must be full-dimensional. The cost is one LP per face certificate plus one per
generator batch, so it scales with the size of the *projected* cone, not the source cone.
"""
function projectcone(P::AbstractMatrix{<:LrsNum}, A::AbstractMatrix{<:LrsNum};
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
    # arithmetic — a float oracle can dispute within tolerance yet return a point that does not
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
    # it grows — probing the remaining facets of a stale hull is counterproductive, their junk
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
    # facets always come from the exact lrs hull: float hull candidates (qhull, see
    # `_qhullconefacets`) were tried and are unstable against the support oracle — over a cone,
    # facet validity is discontinuous in the direction, so noisy candidates read as violated and
    # flood the hull with junk generators, however the tolerance is arranged.
    while true
        F, _, _ = facetsof(_tomatrix(gens, d); rays=true)
        sweep!(certified, F, _normdir) && continue
        rays, _ = extremerays(F, zeros(Rational{Int}, size(F, 1)))
        return F, rays
    end
end

"""
    removeredundancy(A, b; linearity=Int[])

Drop the redundant inequalities from the H-representation `A*x <= b`.

  - returns `(A, b, kept)`, where `kept` indexes the surviving rows of the input
"""
function removeredundancy(A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum}; linearity=Int[])
    normals = Vector{Rational{BigInt}}[]
    bs = Rational{BigInt}[]
    keep = trues(size(A, 1))

    # lrs reports `* n redundant row(s) found:` after the representation, then lists the indices on
    # the next line; that index line is only meaningful once the report has been seen.
    reported = false
    _streamlrs(_ine(A, b; linearity, options=["redund 0 0"]);
               onrow=function (r)
                   push!(normals, -r[2:end])
                   push!(bs, r[1])
               end,
               ontrailer=function (s)
                   if occursin("redundant row", s)
                       reported = true
                   elseif reported && _isnumericrow(s)
                       for t in eachsplit(s)
                           keep[parse(Int, t)] = false
                       end
                       reported = false
                   end
               end)

    return _tomatrix(normals, size(A, 2)), bs, findall(keep)
end
