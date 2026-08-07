const LrsNum = Union{Integer,Rational}

function _writerows(io, rows)
    for r in rows
        for x in r
            print(io, ' ')
            x isa Rational ? print(io, numerator(x), '/', denominator(x)) : print(io, x)
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

# lrs V-representation: each row is [flag v...], flag 0 for a ray and 1 for a vertex. `rays` is
# either a single flag for every row or one per row, so unbounded polyhedra can be described too.
function _ext(V::AbstractMatrix{<:LrsNum}; rays=true, options=String[])
    flags = rays isa Bool ? fill(rays, size(V, 1)) : collect(Bool, rays)
    length(flags) == size(V, 1) ||
        throw(DimensionMismatch("`rays` has $(length(flags)) entries but `V` has $(size(V, 1)) rows"))
    io = IOBuffer()
    println(io, "crafts\nV-representation\nbegin")
    println(io, size(V, 1), " ", size(V, 2) + 1, " rational")
    _writerows(io, (vcat(flags[i] ? 0 : 1, V[i, :]) for i in axes(V, 1)))
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
    open(pipeline(IOBuffer(input), `$(lrs())`)) do io
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
    extremerays(A, b; linearity=Int[])

Convert the H-representation `A*x <= b` into its extreme rays and vertices.

  - `linearity`: indices of rows of `A` that hold with equality
  - returns `(rays, vertices)`, both matrices with one generator per row
"""
function extremerays(A::AbstractMatrix{<:LrsNum}, b::AbstractVector{<:LrsNum}; linearity=Int[])
    rays = Vector{Rational{BigInt}}[]
    verts = Vector{Rational{BigInt}}[]
    _streamlrs(_ine(A, b; linearity);
               onrow=r -> push!(iszero(r[1]) ? rays : verts, r[2:end]))
    d = size(A, 2)
    return _tomatrix(rays, d), _tomatrix(verts, d)
end

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

"""
    facetsof(V; rays=true)

Convert the V-representation `V` into facets, with the incidence of generators on each facet.

  - `V`: one generator per row
  - `rays`: whether the rows are rays of a cone rather than vertices of a polytope
  - returns `(A, b, incidence)` with `A*x <= b`, and `incidence[i, j]` true when generator `j` lies
    on facet `i`
"""
function facetsof(V::AbstractMatrix{<:LrsNum}; rays=true)
    normals = Vector{Rational{BigInt}}[]
    b = Rational{BigInt}[]
    _streamlrs(_ext(V; rays); onrow=function (r)
                   push!(normals, -r[2:end])
                   push!(b, r[1])
               end)
    A = _tomatrix(normals, size(V, 2))
    return A, b, _incidence(A, b, V)
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
